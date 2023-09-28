///////////////////////////////////////////////////////////////////////////////
//
// File: MetricExecutionTime.cpp
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
// Description: Implementation of the execution time metric. A test will fail
// if the execution time of the test falls outside of an assigned tolerance.
//
///////////////////////////////////////////////////////////////////////////////

#include <MetricExecutionTime.h>

#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace Nektar
{
string MetricExecutionTime::type = GetMetricFactory().RegisterCreatorFunction(
    "EXECUTIONTIME", MetricExecutionTime::create);

/**
 * @brief Constructor.
 */
MetricExecutionTime::MetricExecutionTime(TiXmlElement *metric, bool generate)
    : Metric(metric, generate)
{
    // Set up the regular expression to match against.
    m_regex = "^.*Total Computation Time\\s*=\\s*(\\d+\\.?\\d*)s";
    // Inform the Tester this metric supports averaging data from multiple runs.
    m_average = true;
    m_match   = MetricExecutionTimeFieldValue("nil");

    // Parse matching values if not generating.
    if (m_generate)
    {
        return;
    }

    TiXmlElement *value = metric->FirstChildElement("value");
    ASSERTL0(value || m_generate, "Missing value tag for metric.");

    bool hostnameMatch = false;
    while (value)
    {
        ASSERTL0(value->Attribute("tolerance"),
                 "Missing tolerance in execution time metric.");
        ASSERTL0(value->Attribute("hostname"),
                 "Missing hostname in execution time metric.");
        ASSERTL0(!EmptyString(value->GetText()),
                 "Missing value in execution time metric.");

        // Only use the value if it matches the runner's hostname.
        if (value->Attribute("hostname") == boost::asio::ip::host_name())
        {
            MetricExecutionTimeFieldValue val;
            val.m_value     = value->GetText();
            val.m_tolerance = atof(value->Attribute("tolerance"));
            m_match         = val;

            hostnameMatch = true;

            // No use finding any more values
            break;
        }

        value = value->NextSiblingElement("value");
    }

    // If no hostname match was found, pass a skip flag to the test.
    if (!hostnameMatch)
    {
        cerr << "WARNING: No execution time value provided for host "
             << boost::asio::ip::host_name() << ". Skipping metric." << endl;
        m_match.m_skip = true;
    }
}

/**
 * @brief Test output against a regular expression and its expected value.
 */
bool MetricExecutionTime::v_Test(istream &pStdout, istream &pStderr)
{
    boost::ignore_unused(pStdout, pStderr);

    bool success = true;

    // Select istream to use.
    istream &is = m_useStderr ? pStderr : pStdout;

    boost::cmatch matches;
    string line;
    // Vector of execution times found in the output.
    vector<double> times;
    bool matched = false;

    // Process output file line by line searching for regex matches
    while (getline(is, line))
    {

        // Test to see if we have a match on this line.
        if (boost::regex_match(line.c_str(), matches, m_regex))
        {
            // If no matches are found then throw an error.
            if (matches.size() == 1)
            {
                cerr << "No test sections in regex!" << endl;
                return false;
            }

            // Check each regex capture group in turn
            for (int i = 1; i < matches.size(); ++i)
            {
                // Extract the captured string
                string match(matches[i].first, matches[i].second);
                double val;
                try
                {
                    val = boost::lexical_cast<double>(match);

                    // Add value to the list of found execution times.
                    times.push_back(val);
                    matched = true;
                }
                catch (boost::bad_lexical_cast &e)
                {
                    cerr << "Could not convert match " << match << " to double"
                         << endl;
                    success = false;
                    continue;
                }
            }
        }
    }

    if (!matched)
    {
        cerr << "No execution times were found in the test output." << endl;
        success = false;
    }

    // Average the found execution times.
    double avgTime = 0.0;
    for (unsigned int i = 0; i < times.size(); ++i)
    {
        avgTime += times[i];
    }
    avgTime /= times.size();

    // If a skip flag is set (due to no matching hostname) then return
    // successful test and output the average time.
    if (m_match.m_skip)
    {
        cerr << "Average execution time for host "
             << boost::asio::ip::host_name() << ": " << avgTime << endl;
        return success;
    }

    ASSERTL0(m_match.m_value != "nil",
             "No test conditions defined for execution time.");

    // Check that the average time is within the tolerance.
    if (fabs(avgTime - boost::lexical_cast<double>(m_match.m_value)) >
            m_match.m_tolerance &&
        !m_match.m_skip)
    {
        cerr << endl;
        cerr << "Failed tolerance match." << endl;
        cerr << "  Expected avg execution time: " << m_match.m_value << " +/- "
             << m_match.m_tolerance << endl;
        cerr << "  Actual avg: " << avgTime << endl;

        for (unsigned int i = 0; i < times.size(); ++i)
        {
            cerr << "       Run " << i << ": " << times[i] << endl;
        }

        success = false;
    }

    return success;
}

/**
 * @brief Generate accepted values for this metric using an output or error
 * file.
 */
void MetricExecutionTime::v_Generate(istream &pStdout, istream &pStderr)
{
    boost::ignore_unused(pStderr);

    // Select istream to use
    istream &is = m_useStderr ? pStderr : pStdout;

    boost::cmatch matches;

    string line;
    // Vector of execution times found in the output
    vector<double> sampleTimes;
    bool matched = false;

    // Process output file line by line searching for regex matches
    while (getline(is, line))
    {
        // Test to see if we have a match on this line
        if (boost::regex_match(line.c_str(), matches, m_regex))
        {
            // If no fields in regex then throw an error
            ASSERTL0(matches.size() != 1, "No test sections in regex!");

            // Check each regex capture group in turn
            for (int i = 1; i < matches.size(); ++i)
            {
                // Extract the captured string
                string match(matches[i].first, matches[i].second);
                double val;
                try
                {
                    val = boost::lexical_cast<double>(match);

                    sampleTimes.push_back(val);
                    matched = true;
                }
                catch (boost::bad_lexical_cast &e)
                {
                    cerr << "Could not convert match " << match << " to double"
                         << endl;
                    continue;
                }
            }
        }
    }

    // Stores the accepted average execution time and a default tolerance.
    MetricExecutionTimeFieldValue okValue;

    if (matched)
    {
        // Average the found execution times and set it as the accepted value.
        double avgTime = 0.0;
        for (unsigned int i = 0; i < sampleTimes.size(); ++i)
        {
            avgTime += sampleTimes[i];
        }
        avgTime /= sampleTimes.size();

        okValue.m_value = to_string(avgTime);
        m_match         = okValue;
    }
    else
    {
        cerr << "No execution times were found in the test output. Value set "
                "to 'nil'."
             << endl;
    }

    // If we are not a derived class then create a new structure.
    if (m_type == "EXECUTIONTIME")
    {
        // Remove values if they already exist.
        while (m_metric->FirstChildElement("value"))
        {
            ASSERTL0(
                m_metric->RemoveChild(m_metric->FirstChildElement("value")),
                "Couldn't remove value from metric!");
        }

        TiXmlElement *val = new TiXmlElement("value");

        val->SetAttribute("tolerance",
                          boost::lexical_cast<string>(m_match.m_tolerance));
        val->SetAttribute("hostname", boost::asio::ip::host_name());
        val->LinkEndChild(new TiXmlText(m_match.m_value));
        m_metric->LinkEndChild(val);
    }
}
} // namespace Nektar
