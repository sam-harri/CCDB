///////////////////////////////////////////////////////////////////////////////
//
// File: ThreadStd.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#define NOMINMAX

#include <algorithm>
#include <iostream>

#include <LibUtilities/BasicUtils/ThreadStd.h>

namespace Nektar::Thread
{

std::string ThreadManagerStd::className =
    GetThreadManagerFactory().RegisterCreatorFunction(
        "ThreadManagerStd", ThreadManagerStd::Create,
        "Threading using std::thread.");

/**
 * @param numWorkers The number of threads to start (including master thread).
 * @note Do not use, use factory instead.
 */
ThreadManagerStd::ThreadManagerStd(unsigned int numT)
    : m_numThreads(numT), m_numWorkers(numT), m_masterQueue(),
      m_masterQueueMutex(), m_masterActiveMutex(), m_masterQueueCondVar(),
      m_masterActiveCondVar(), m_chunkSize(1), m_schedType(e_dynamic),
      m_threadMap()
{
    using namespace std;
    try
    {
        m_threadList       = new ThreadWorkerStd *[m_numThreads];
        m_threadThreadList = new std::thread *[m_numThreads];
        m_threadBusyList   = new bool[m_numThreads];
        m_threadActiveList = new bool[m_numThreads];
    }
    catch (exception &e)
    {
        cerr << "Exception while allocating thread storage: " << e.what()
             << endl;
        abort();
    }
    unsigned int i = 0;
    while (i < m_numThreads)
    {
        ThreadWorkerStd *tw;
        try
        {
            tw = new ThreadWorkerStd(this, i);
        }
        catch (exception &e)
        {
            cerr << "Exception while allocating worker threads: " << e.what()
                 << endl;
            abort();
        }

        m_threadList[i]       = tw;
        m_threadBusyList[i]   = false;
        m_threadActiveList[i] = true;

        try
        {
            m_threadThreadList[i] = new std::thread(std::ref(*tw));
            std::thread::id id    = m_threadThreadList[i]->get_id();
            m_threadMap[id]       = i;
        }
        catch (...)
        {
            std::cerr << "Exception while creating worker threads" << std::endl;
            abort();
        }
        i++;
    }
    m_masterThreadId = std::this_thread::get_id();
    m_barrier        = new Barrier(m_numWorkers > 0 ? m_numWorkers : 1);
    m_type           = "Threading with Std";
}

/**
 * Terminates all running threads (they will finish their current job),
 * releases resources and destructs.
 */
ThreadManagerStd::~ThreadManagerStd()
{
    // This is an immediate teardown.  We attempt to kill everything.
    // we daren't lock anything as we may cause a deadlock
    for (unsigned int i = 0; i < m_numThreads; i++)
    {
        m_threadList[i]->Stop();
    }

    m_masterQueueCondVar.notify_all();
    m_masterActiveCondVar.notify_all();
    for (unsigned int i = 0; i < m_numThreads; i++)
    {
        m_threadThreadList[i]->join();
        delete m_threadThreadList[i];
        delete m_threadList[i];
    }

    delete[] m_threadList;
    delete[] m_threadThreadList;
    delete[] m_threadActiveList;
    delete[] m_threadBusyList;
    delete m_barrier;
}

/**
 *
 */
void ThreadManagerStd::v_QueueJobs(std::vector<ThreadJob *> &joblist)
{
    std::vector<ThreadJob *>::iterator it;
    for (it = joblist.begin(); it < joblist.end(); ++it)
    {
        QueueJob(*it);
    }
}

/*
 *
 */
void ThreadManagerStd::v_QueueJob(ThreadJob *job)
{
    Lock masterQueueLock(m_masterQueueMutex); // locks the queue
    m_masterQueue.push(job);
    m_masterQueueCondVar.notify_all(); // alert a waiting thread.
} // queue unlocked

/**
 *
 */
bool ThreadManagerStd::IsWorking()
{
    bool working = false;
    Lock masterActiveLock(m_masterActiveMutex);
    for (unsigned int i = 0; i < m_numWorkers; i++)
    {
        working = working || m_threadBusyList[i];
    }
    return working;
}

/**
 *
 */
void ThreadManagerStd::v_SetChunkSize(unsigned int chnk)
{
    Lock masterQueueLock(m_masterQueueMutex); // locks the queue
    m_chunkSize = std::max(chnk, 1U);
}

/**
 *
 */
void ThreadManagerStd::v_SetSchedType(SchedType s)
{
    Lock masterQueueLock(m_masterQueueMutex); // locks the queue
    m_schedType = s;
}

/**
 *
 */
bool ThreadManagerStd::v_InThread()
{
    std::thread::id id = std::this_thread::get_id();
    return (id != m_masterThreadId);
}

/**
 *
 */
void ThreadManagerStd::v_Wait()
{
    bool working;
    Lock masterQueueLock(m_masterQueueMutex); // locks the queue
    working = IsWorking();
    m_masterActiveCondVar.notify_all();
    m_masterQueueCondVar.notify_all();
    while (!m_masterQueue.empty() || working)
    {
        // while waiting, master queue is unlocked
        m_masterQueueCondVar.wait(masterQueueLock);
        // on exiting wait master queue is locked again
        working = IsWorking();
    }
}

/**
 *
 */
unsigned int ThreadManagerStd::v_GetNumWorkers()
{
    return m_numWorkers;
}

/**
 *
 */
unsigned int ThreadManagerStd::v_GetWorkerNum()
{
    std::thread::id id = std::this_thread::get_id();
    return m_threadMap[id];
}

/**
 *
 */
void ThreadManagerStd::SetNumWorkersImpl(const unsigned int num)
{
    Lock masterActiveLock(m_masterActiveMutex); // locks the active

    if (m_numWorkers == num)
    {
        return;
    }

    delete m_barrier;
    m_barrier = new Barrier(num > 0 ? num : 1);

    m_numWorkers = num;
    for (unsigned int i = 0; i < m_numThreads; i++)
    {
        m_threadActiveList[i] = i < m_numWorkers ? true : false;
    }
    m_masterActiveCondVar.notify_all();
} // Lock on active released here

/**
 *
 */
void ThreadManagerStd::v_SetNumWorkers(const unsigned int num)
{
    unsigned int n;
    n = std::min(num, m_numThreads);
    n = std::max(n, static_cast<unsigned int>(0));
    SetNumWorkersImpl(n);
}

/**
 *
 */
void ThreadManagerStd::v_SetNumWorkers()
{
    SetNumWorkersImpl(m_numThreads);
}

/**
 *
 */
unsigned int ThreadManagerStd::v_GetMaxNumWorkers()
{
    return m_numThreads;
}

/**
 *
 */
void ThreadManagerStd::v_Hold()
{
    m_barrier->Wait();
}

/**
 *
 */
const std::string &ThreadManagerStd::v_GetType() const
{
    return m_type;
}

/**
 * @param threadManager Pointer to the ThreadManagerStd that is controlling
 *                      this worker.
 * @param workerNum Unique number from 0..(number_of_threads - 1)
 *
 * Called by the ThreadManagerStd instance.
 */
ThreadWorkerStd::ThreadWorkerStd(ThreadManagerStd *tm, unsigned int workerNum)
    : m_threadManager(tm), m_workerQueue(), m_keepgoing(true),
      m_threadNum(workerNum)
{
    // Nothing to see here
}

/**
 * Winds up this thread's execution.  Jobs in its queue are lost.
 */
ThreadWorkerStd::~ThreadWorkerStd()
{
    if (m_keepgoing)
    {
        std::cerr << "Warning: ThreadWorker: " << m_threadNum
                  << "destroyed while running!" << std::endl;
    }
    // on destuction the m_workerQueue will be destructed and that
    // will destruct any ThreadJobs still in there.
}

/**
 *
 */
void ThreadWorkerStd::LoadJobs()
{
    // Lock the master queue
    Lock masterQueueLock(m_threadManager->m_masterQueueMutex);
    m_threadManager->m_threadBusyList[m_threadNum] = false;
    m_threadManager->m_masterQueueCondVar.notify_all();
    while (m_threadManager->m_masterQueue.empty() && m_keepgoing)
    {
        // while waiting, master queue is unlocked
        m_threadManager->m_masterQueueCondVar.wait(masterQueueLock);
        // on exiting wait master queue is locked again
    }
    bool active;
    {
        Lock masterActiveLock(m_threadManager->m_masterActiveMutex);
        active = m_threadManager->m_threadActiveList[m_threadNum];
    }
    if (active && m_keepgoing)
    {
        unsigned int numToLoad = GetNumToLoad();
        while (m_workerQueue.size() < numToLoad &&
               !m_threadManager->m_masterQueue.empty())
        {
            ThreadJob *tj = m_threadManager->m_masterQueue.front();
            m_workerQueue.push(tj);
            m_threadManager->m_masterQueue.pop();
        }
    }
    m_threadManager->m_threadBusyList[m_threadNum] = !m_workerQueue.empty();
} // lock on master queue released here

/**
 *
 */
unsigned int ThreadWorkerStd::GetNumToLoad()
{
    unsigned int numToLoad = 0;
    switch (m_threadManager->m_schedType)
    {
        case e_guided:
            numToLoad = std::max(
                static_cast<unsigned long>(m_threadManager->m_chunkSize),
                static_cast<unsigned long>(
                    m_threadManager->m_masterQueue.size() /
                    (2 * m_threadManager->m_numWorkers + 1)));
            break;

        case e_dynamic:
            numToLoad = m_threadManager->m_chunkSize;
            break;

        default:
            NEKERROR(ErrorUtil::efatal, "Invalid value for SchedType.");
            break;
    }
    return numToLoad;
}

/**
 *
 */
void ThreadWorkerStd::WaitForActive()
{
    Lock masterActiveLock(m_threadManager->m_masterActiveMutex);

    while (!m_threadManager->m_threadActiveList[m_threadNum] && m_keepgoing)
    {
        // while waiting, master active is unlocked
        m_threadManager->m_masterActiveCondVar.wait(masterActiveLock);
        // on exiting wait master active is locked again
    }
}

/**
 *
 */
void ThreadWorkerStd::MainLoop()
{
    while (m_keepgoing)
    {
        WaitForActive();
        LoadJobs();
        RunJobs();
    }
} // exiting here should terminate the thread

/**
 *
 */
void ThreadWorkerStd::RunJobs()
{
    while (!m_workerQueue.empty() && m_keepgoing)
    {
        ThreadJob *tj;
        try
        {
            tj = m_workerQueue.front();
            tj->SetWorkerNum(m_threadNum);
            tj->Run();
            m_workerQueue.pop();
            delete tj;
        }
        catch (...)
        {
            // something bad happened, probably time to die
            // maybe signal ThreadManager
            throw;
        }
    }
}

} // namespace Nektar::Thread
