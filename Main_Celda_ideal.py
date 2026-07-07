"""
Main code where the principle flow is controlled
"""

import os
import subprocess
import sys

file_name = "CeldillaV1_ideal"  # cell to study
reports_dir = r"D:\Eduardo_Rodriguez\Python_archivos\Bucle\Reportes"  # directory to save all results
temporal = r"D:\Eduardo_Rodriguez\Temporal" # temporary directory to store and delete obsolete files

# Dictionary with every code and its number
CODE_MAP = {
    "1": "Calculo_convergencia_uniaxial.py",
    "2": "Calculo_cortante_elastico.py",
    "3": "Calculo_uniaxial_plastico.py",
    "4": "Calculo_cortante_plastico.py",
    "5": "Criterio_fallo.py",
}

# --- User menu ---
print("--- Simulation Selector to execute ---")
print("1 - Calculo_convergencia_uniaxial")
print("2 - Calculo_cortante_elastico")
print("3 - Calculo_uniaxial_plastico")
print("4 - Calculo_cortante_plastico")
print("5 - Criterio_fallo")
print("-" * 45)

selection = input(
    "Enter the numbers of the codes you want to run (e.g., '12' or '135'): "
).strip()

print(f"\nStarting execution sequence for model: {file_name}...\n")

# Loop of simulations chosen
for digit in selection:
    if digit in CODE_MAP:
        script_to_run = CODE_MAP[digit]
        print(f"==================================================")
        print(f"Executing [{digit}]: {script_to_run}...")
        print(f"==================================================")

        try:
            # Executes the script sending 'file_name' and 'reports_dir' as arguments
            result = subprocess.run(
                [sys.executable, script_to_run, file_name, reports_dir, temporal],
                check=True,  # Stops the main script if the sub-code fails
                text=True,  # Captures console output as text
            )
            print(f"\n[OK] {script_to_run} finished successfully.\n")

        except subprocess.CalledProcessError:
            print(
                f"\n[ERROR] The script {script_to_run} failed. Stopping the sequence."
            )
            break
    else:
        print(
            f"Warning: The digit '{digit}' does not correspond to any valid simulation. Skipping..."
        )

print("Sequential process finished.")