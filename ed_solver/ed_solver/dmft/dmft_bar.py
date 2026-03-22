RED     = "\033[31m"
GREEN   = "\033[32m"
YELLOW  = "\033[33m"
RESET   = "\033[0m"

def bar_fmt(curr, prev):
    if curr <= prev:
        return f"{GREEN}{curr:.2e}{RESET}"
    else:
        return f"{RED}{curr:.2e}{RESET}"\
        