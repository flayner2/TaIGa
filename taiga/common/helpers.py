import sys
import logging as log


def sanitize_version():
    if sys.version_info < (3, 6):
        sys.exit("""You are running TaIGa with a Python version older than 3.6x.\
        Please, try again with 'python3 TaIGa.py' or 'python3.6 TaIGa.py'.""")


def sanitize_args(args):
    if ((args.multi and (args.same or args.single or args.tid))
       or (args.same and (args.multi or args.single or args.tid))
       or (args.single and (args.same or args.multi or args.tid))
       or (args.tid and (args.same or args.single or args.multi))):
        log.error("\nERROR: Please run TaIGa with only one of the possible input type optional arguments.")
        sys.exit()


def config_log(verbose):
    if verbose:
        log.basicConfig(
            filename='TaIGa_run.log', format="%(message)s", level=log.DEBUG)
    else:
        log.basicConfig(format="%(message)s", level=log.DEBUG)


def sanitize_name_list(correction, tid, missing_corrected, missing_taxid, input_names):
    if (correction) and (not tid):
        # Check for the organisms with missing corrected name or taxid and remove their name form the names list
        for missed in missing_corrected:
            if missed in input_names:
                input_names.remove(missed)
    if not tid:
        for missed in missing_taxid:
            if missed in input_names:
                input_names.remove(missed)
