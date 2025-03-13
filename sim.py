# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "jinja2",
# ]
# ///

import logging
import os
import glob
import shutil
import subprocess
import json
from pathlib import Path
from jinja2 import Environment, FileSystemLoader


# Configuration
class Config:
    # path params
    NAME = "TEST"
    SESSION_TEMPLATE = "session_template.xml"
    MOCK_SIM = "../../mock_sim.sh"

    # Phase 1 params
    PHASE1_DT = 6.25e-9
    PHASE1_STEPS = 10

    # Phase 2 params
    PHASE2_INITIAL_DT = 7e-7
    PHASE2_DESIRED_RUNTIME = 6e-5
    PHASE2_CFL_TARGET = 0.2
    PHASE2_CFL_MIN = 0.195
    PHASE2_CFL_MAX = 0.205

    # Phase 3 params
    PHASE3_DESIRED_RUNTIME = 6e-4
    # Simulation params
    VELOCITY_NUM_MODES = 3
    PRESSURE_NUM_MODES = 2


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


class SimulationManager:
    def __init__(self):
        self.base_dir = Path(Config.NAME)
        self.phase_dirs = {
            1: self.base_dir / "phase1",
            2: self.base_dir / "phase2",
            3: self.base_dir / "phase3",
        }
        self.env = Environment(loader=FileSystemLoader("."))
        self.current_dt = 0
        logger.info(
            f"Initialized simulation manager with base directory: {self.base_dir}"
        )

    def setup_directories(self):
        """Create necessary directories for the simulation."""
        for phase, phase_dir in self.phase_dirs.items():
            phase_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created/verified directory for phase {phase}: {phase_dir}")

    def find_latest_checkpoint(self, phase_dir):
        """Find the latest checkpoint file in the given directory."""
        checkpoint_files = glob.glob(str(phase_dir / "coil_*.chk"))
        if not checkpoint_files:
            logger.warning(f"No checkpoint files found in {phase_dir}")
            return None
        latest = max(
            checkpoint_files, key=lambda x: int(x.split("_")[-1].split(".")[0])
        )
        logger.info(f"Found latest checkpoint in {phase_dir}: {latest}")
        return latest

    def copy_checkpoint_as_restart(self, source_dir, target_dir):
        """Copy the latest checkpoint as a restart file to the target directory."""
        latest_checkpoint = self.find_latest_checkpoint(source_dir)
        if latest_checkpoint:
            restart_file = target_dir / "restart.rst"
            shutil.copy2(latest_checkpoint, restart_file)
            logger.info(
                f"Copied checkpoint {latest_checkpoint} to restart file {restart_file}"
            )
            return restart_file.name
        logger.warning(f"No checkpoint found to copy from {source_dir} to {target_dir}")
        return None

    def update_session_file(self, phase, dt, runtime, num_steps, initial_conditions_file=None):
        """Update the session file with the given parameters."""
        template = self.env.get_template(Config.SESSION_TEMPLATE)
        content = template.render(
            timestep=dt,
            velocity_num_modes=Config.VELOCITY_NUM_MODES,
            pressure_num_modes=Config.PRESSURE_NUM_MODES,
            runtime=runtime,
            num_steps=num_steps,
            initial_conditions_file=initial_conditions_file,
        )

        session_file = self.phase_dirs[phase] / "session.xml"
        with open(session_file, "w") as f:
            f.write(content)
        logger.info(f"Updated session file for phase {phase}: {session_file}")
        if initial_conditions_file:
            logger.info(
                f"Session file configured to use initial conditions from: {initial_conditions_file}"
            )

    def run_phase_one(self):
        """Run phase 1 simulation with small timestep."""
        logger.info(f"Starting phase 1 in directory: {self.phase_dirs[1]}")
        logger.info(
            f"Running phase 1 with dt={Config.PHASE1_DT}, steps={Config.PHASE1_STEPS}"
        )
        cmd = [Config.MOCK_SIM, str(Config.PHASE1_DT), str(Config.PHASE1_STEPS)]
        logger.info(f"Command: {' '.join(cmd)}")

        # Run the simulation and capture output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=self.phase_dirs[1],
        )

        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line and line[0] == "{":
                data = json.loads(line)
                cfl = data["CFL"]
                logger.info(f"Phase 1 CFL: {cfl}")

        process.wait()
        if process.returncode != 0:
            logger.error(
                f"Phase 1 simulation failed with return code {process.returncode}"
            )
            raise RuntimeError("Phase 1 simulation failed")
        logger.info("Phase 1 completed successfully")

    def run_phase_two(self):
        """Run phase 2 simulation with adaptive timestep based on CFL."""
        restart_file = self.copy_checkpoint_as_restart(
            self.phase_dirs[1], self.phase_dirs[2]
        )
        logger.info(f"Starting phase 2 in directory: {self.phase_dirs[2]}")
        self.current_dt = Config.PHASE2_INITIAL_DT
        iteration = 0

        while True:  # Keep trying until we get the right CFL
            iteration += 1
            # Calculate rounded dt and optimal steps
            rounded_dt, num_steps, actual_runtime = self.calculate_steps(
                self.current_dt, Config.PHASE2_DESIRED_RUNTIME
            )
            
            self.update_session_file(
                phase=2,
                dt=rounded_dt,
                runtime=actual_runtime,
                num_steps=num_steps,
                initial_conditions_file=restart_file,
            )
            
            logger.info(f"Phase 2 iteration {iteration} with dt={rounded_dt:.2e}")
            cmd = [Config.MOCK_SIM, str(rounded_dt), str(num_steps)]
            logger.info(f"Command: {' '.join(cmd)}")

            # Run the simulation and capture output
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=self.phase_dirs[2],
            )

            line = None
            while True:
                line = process.stdout.readline()
                if line and line[0] == "{":
                    break

            data = json.loads(line)
            cfl = data["CFL"]
            logger.info(
                f"Phase 2 CFL of {cfl} at iteration {iteration} with dt={rounded_dt:.2e}"
            )

            # If CFL is out of range, adjust timestep and restart
            if cfl < Config.PHASE2_CFL_MIN or cfl > Config.PHASE2_CFL_MAX:
                old_dt = self.current_dt
                self.current_dt = self.adjust_timestep(self.current_dt, cfl)
                logger.info(
                    f"CFL {cfl} out of range [{Config.PHASE2_CFL_MIN}, {Config.PHASE2_CFL_MAX}]. "
                    f"Adjusting timestep from {old_dt:.2e} to {self.current_dt:.2e} and restarting"
                )
                process.terminate()
                process.wait()
                continue

            # If we get here, CFL is good, let the simulation complete
            # remove all the bad checkpoints
            for file in glob.glob(str(self.phase_dirs[2] / "coil_*.chk")):
                os.remove(file)
            process.wait()
            if process.returncode != 0:
                logger.error(
                    f"Phase 2 simulation failed with return code {process.returncode}"
                )
                raise RuntimeError("Phase 2 simulation failed")

            logger.info(f"Phase 2 completed successfully with dt={rounded_dt:.2e}")
            break  # Exit the loop once we complete a successful run

    def run_phase_three(self):
        """Run phase 3 simulation with final timestep."""
        logger.info(f"Starting phase 3 in directory: {self.phase_dirs[3]}")
        
        # Calculate rounded dt and optimal steps
        rounded_dt, num_steps, actual_runtime = self.calculate_steps(
            self.current_dt, Config.PHASE3_DESIRED_RUNTIME
        )
        
        logger.info(f"Running phase 3 with dt={rounded_dt:.2e}, steps={num_steps}")
        cmd = [Config.MOCK_SIM, str(rounded_dt), str(num_steps)]
        logger.info(f"Command: {' '.join(cmd)}")

        # Run the simulation and capture output
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=self.phase_dirs[3],
        )

        # Process output in real-time
        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line and line[0] == "{":
                data = json.loads(line)
                logger.info(f"Phase 3 CFL: {data['CFL']}")

        process.wait()
        if process.returncode != 0:
            logger.error(
                f"Phase 3 simulation failed with return code {process.returncode}"
            )
            raise RuntimeError("Phase 3 simulation failed")
        logger.info("Phase 3 completed successfully")

    def adjust_timestep(self, current_dt, current_cfl):
        """Adjust timestep based on CFL number using law of similitude."""
        if current_cfl > Config.PHASE2_CFL_TARGET:
            return current_dt * (Config.PHASE2_CFL_TARGET / current_cfl)
        else:
            return current_dt * (Config.PHASE2_CFL_TARGET / current_cfl)

    def round_timestep(self, dt):
        """Round timestep to 2 decimal places in scientific notation."""
        # Convert to scientific notation with 2 decimal places
        return float(f"{dt:.2e}")

    def calculate_steps(self, dt, desired_runtime):
        """Calculate number of steps that gives closest runtime to desired."""
        # Round the timestep first
        rounded_dt = self.round_timestep(dt)
        
        # Calculate exact number of steps needed
        exact_steps = desired_runtime / rounded_dt
        
        # Round to nearest integer
        num_steps = round(exact_steps)
        
        # Calculate actual runtime with this number of steps
        actual_runtime = rounded_dt * num_steps
        
        logger.info(f"Timestep rounded from {dt:.2e} to {rounded_dt:.2e}")
        logger.info(f"Desired runtime: {desired_runtime:.2e}, Actual runtime: {actual_runtime:.2e}")
        logger.info(f"Using {num_steps} steps")
        
        return rounded_dt, num_steps, actual_runtime

    def phase_one(self):
        """Run phase 1 with small timestep."""
        logger.info("Starting Phase 1")
        self.setup_directories()
        self.update_session_file(
            1, Config.PHASE1_DT, runtime=Config.PHASE1_DT * Config.PHASE1_STEPS, num_steps=Config.PHASE1_STEPS
        )
        self.run_phase_one()

    def phase_two(self):
        """Run phase 2 with adaptive timestep based on CFL."""
        logger.info("Starting Phase 2")

        self.run_phase_two()

    def phase_three(self):
        """Run phase 3 with final timestep from phase 2."""
        logger.info("Starting Phase 3")

        # Copy final checkpoint from phase 2 as restart
        restart_file = self.copy_checkpoint_as_restart(
            self.phase_dirs[2], self.phase_dirs[3]
        )

        # Calculate rounded dt and optimal steps
        rounded_dt, num_steps, actual_runtime = self.calculate_steps(
            self.current_dt, Config.PHASE3_DESIRED_RUNTIME
        )
        
        logger.info(f"Starting phase 3 with timestep from phase 2: {rounded_dt:.2e}")

        # Create session file with restart configuration
        self.update_session_file(
            phase=3,
            dt=rounded_dt,
            runtime=actual_runtime,
            num_steps=num_steps,
            initial_conditions_file=restart_file,
        )

        self.run_phase_three()

    def run(self):
        """Executes all simulation phases."""
        logger.info(f"Starting simulation in base directory: {self.base_dir}")
        self.phase_one()
        self.phase_two()
        self.phase_three()
        logger.info("All phases completed successfully")


if __name__ == "__main__":
    manager = SimulationManager()
    manager.run()
