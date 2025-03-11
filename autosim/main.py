import os
import subprocess
import time
import shutil
import re
import numpy as np
from jinja2 import Environment, FileSystemLoader

# Configuration
SESSION_TEMPLATE = "session_template.xml"
SESSION_FILE = "session.xml"
CHECKPOINT_DIR = "checkpoint/"
NUM_PROCESSES = 4
MPI_CMD = "mpiexec -np {num_processes} ../build/solvers/IncNavierStokesSolver/IncNavierStokesSolver coil.xml session.xml -v"
TARGET_CFL = 0.2
INITIAL_DT = 6.25e-9
INITIAL_CFL = 0.01

class SimulationManager:
    def __init__(self):
        self.current_cfl = INITIAL_CFL
        self.current_dt = INITIAL_DT
        self.stage = 1
        self.running_process = None

    def update_session_file(self, dt=None, restart_file=None):
        """Updates the session file using Jinja2 template."""
        env = Environment(loader=FileSystemLoader('.'))
        template = env.get_template(SESSION_TEMPLATE)
        rendered = template.render(timestep=dt, restart_file=restart_file, stage=self.stage)
        with open(SESSION_FILE, "w") as f:
            f.write(rendered)

    def get_latest_checkpoint(self):
        """Returns the latest checkpoint file."""
        files = [f for f in os.listdir(CHECKPOINT_DIR) if f.endswith(".chk")]
        if not files:
            return None
        latest_chk = max(files, key=lambda f: int(f.split("_")[-1].split(".")[0]))
        return os.path.join(CHECKPOINT_DIR, latest_chk)

    def extract_cfl_from_output(self, output_line):
        """Extracts CFL value from a line of standard output."""
        match = re.search(r"CFL:\s*([\d\.]+)", output_line)
        return float(match.group(1)) if match else None

    def adjust_timestep(self, last_cfl):
        """Adjusts the timestep using similarity principle."""
        new_dt = (TARGET_CFL / last_cfl) * self.current_dt
        print(f"[INFO] Adjusting timestep: Old={self.current_dt:.6e}, New={new_dt:.6e}, CFL={last_cfl:.4f}")
        self.current_dt = new_dt
        self.restart_simulation()

    def restart_simulation(self):
        """Stops the running process, updates the session file, and restarts the simulation."""
        if self.running_process:
            self.running_process.terminate()
            self.running_process.wait()

        latest_chk = self.get_latest_checkpoint()
        if latest_chk:
            rst_file = latest_chk.replace(".chk", ".rst")
            shutil.copy(latest_chk, rst_file)
            self.update_session_file(dt=self.current_dt, restart_file=rst_file)

        self.run_simulation(real_time_adjust=True)

    def run_simulation(self, real_time_adjust=False):
        """Runs the simulation and monitors output for CFL adjustments."""
        cmd = MPI_CMD.format(num_processes=NUM_PROCESSES)
        print(f"[INFO] Starting simulation (Stage {self.stage})...")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        self.running_process = process

        if real_time_adjust:
            for line in iter(process.stdout.readline, ''):
                print(line.strip())
                cfl_value = self.extract_cfl_from_output(line)
                if cfl_value:
                    self.adjust_timestep(cfl_value)
            
        process.wait()

    def phase_one(self):
        """Runs Phase 1: Initialization at low CFL."""
        print("[INFO] Starting Phase 1 (Initialization)...")
        self.stage = 1
        self.current_dt = INITIAL_DT
        self.update_session_file(dt=self.current_dt)
        self.run_simulation()
        print("[INFO] Phase 1 completed. Moving to Phase 2...")

    def phase_two(self):
        """Runs Phase 2: CFL Adjustment with real-time correction."""
        print("[INFO] Starting Phase 2 (CFL Adjustment)...")
        self.stage = 2
        latest_chk = self.get_latest_checkpoint()
        if latest_chk:
            rst_file = latest_chk.replace(".chk", ".rst")
            shutil.copy(latest_chk, rst_file)
            self.update_session_file(dt=self.current_dt, restart_file=rst_file)
        self.run_simulation(real_time_adjust=True)

    def phase_three(self):
        """Runs Phase 3: Final simulation with adjusted CFL."""
        print("[INFO] Starting Phase 3 (Final Simulation)...")
        self.stage = 3
        latest_chk = self.get_latest_checkpoint()
        if latest_chk:
            rst_file = latest_chk.replace(".chk", ".rst")
            shutil.copy(latest_chk, rst_file)
            self.update_session_file(dt=self.current_dt, restart_file=rst_file)
        self.run_simulation()

    def run(self):
        """Executes all simulation phases."""
        self.phase_one()
        self.phase_two()
        self.phase_three()


if __name__ == "__main__":
    manager = SimulationManager()
    manager.run()
