import sys
import os


def create_dir(*args):
    """
    Creates a directory if not exist.
    :param args: Strings to use in os.path.join
    :return: Path of created file
    """
    dir_path = os.path.join(*args)
    if not os.access(dir_path, os.W_OK) or not os.path.isdir(dir_path):  # Create directory if not exist
        os.mkdir(dir_path)
        print(f"{bcolors.WARNING}Directory created: {dir_path}{bcolors.ENDC}")
    return dir_path


def progressBarForTerminal (iteration, total, prefix ='Progress:', suffix ='', decimals = 1, barLength = 50):
    """
    This function should be called inside of loop, gives the loop's progress.
    :param iteration: It is integer. It is current iteration.
    :param total: It is integer. It is total iteration.
    :param prefix: It is string. It will be placed before progress bar.
    :param suffix: It is string. It will be placed after progress bar.
    :param decimals: It is integer. It is number of decimals in percent complete.
    :param barLength: It is integer. It is character length of bar.
    :return: It is void function. Nothing is returned.
    """
    filledLength = int(round(barLength * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


class bcolors:
    HEADER = '\033[95m\033[1m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
