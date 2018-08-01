import shutil

def add_missing_options_to_nml_file(fname, line_start=None):
    # Some of the flags we may wish to change are missin from the default
    # file so we can't adjust them via this script...add them

    if line_start is None:
        line_start = sum(1 for line in open(fname)) - 1

    f = open(fname, "r")
    contents = f.readlines()
    f.close()

    arg = "   cable_user%GS_SWITCH = 'medlyn'\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%GW_MODEL = .FALSE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%or_evap = .TRUE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    tmp_fname = "tmp.nml"
    f = open(tmp_fname, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    shutil.move(tmp_fname, fname)
