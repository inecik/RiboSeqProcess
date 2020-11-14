# Kemal Inecik

import sys
import subprocess
from time import gmtime, strftime, sleep

try:
    while True:
        sys.stdout.flush()
        output = subprocess.check_output(f"ls -l{sys.argv[1]}", shell=True).decode("utf-8")
        line = output.count('\n')
        print('\033[94m' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '\033[0m')
        sys.stdout.write('\033[96m' + output + '\033[0m')
        sys.stdout.write("\033[F" * (line + 1))
        sleep(5)
except:
    print("\n" * (line + 1) + '\033[0m')