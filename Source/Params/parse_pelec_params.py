#!/usr/bin/env python


# This script parses the list of C++ runtime parameters and writes the
# necessary header files to make them available in PeleC's C++ routines
#
# parameters have the format:
#
#   name  type  default
#
# the first three (name, type, default) are mandatory:
#
#   name: the name of the parameter.  This will be the same name as the
#     variable in C++ unless a pair is specified as (name, cpp_name)
#
#   type: the C++ data type (int, Real, string)
#
#   default: the default value.  If specified as a pair, (a, b), then
#     the first value is the normal default and the second is for
#     debug mode (#ifdef AMREX_DEBUG)
#
# Any line beginning with a "#" is ignored
#
# Commands begin with a "@":
#
#    @namespace: sets the namespace that these will be under (see below)
#      it also gives the C++ class name.
#      if we include the keyword "static" after the name, then the parameters
#      will be defined as static member variables in C++
#
#      e.g. @namespace pelec PeleC static
#
# Note: categories listed in the input file aren't used for code generation
# but are used for the documentation generation
#
#
# For a namespace, name, we write out:
#
#   -- name_params.H  (for pelec, included in PeleC.H):
#      declares the static variables of the PeleC class
#
#   -- name_defaults.H  (for pelec, included in peleC.cpp):
#      sets the defaults of the runtime parameters
#
#   -- name_queries.H  (for pelec, included in PeleC.cpp):
#      does the parmparse query to override the default in C++
#

import argparse
import re
import sys

CWARNING = """
// This file is automatically created by parse_pelec_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh\n
"""

param_include_dir = "param_includes/"


class Param(object):
    """the basic parameter class.  For each parameter, we hold the name,
    type, and default.  For some parameters, we also take a second
    value of the default, for use in debug mode (delimited via
    #ifdef AMREX_DEBUG)

    """

    def __init__(
        self,
        name,
        dtype,
        default,
        cpp_var_name=None,
        namespace=None,
        cpp_class=None,
        static=None,
        debug_default=None,
        ifdef=None,
    ):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.cpp_var_name = cpp_var_name

        self.namespace = namespace
        self.cpp_class = cpp_class

        if static is None:
            self.static = 0
        else:
            self.static = static

        self.debug_default = debug_default

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

    def get_default_string(self):
        # this is the line that goes into pelec_defaults.H included
        # into PeleC.cpp

        if self.dtype == "int":
            tstr = "int {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "amrex::Real {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "std::string {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "bool {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "dim_array":
            tstr = "amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> {}::{}".format(
                self.cpp_class, self.cpp_var_name
            )
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if self.ifdef is not None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        if self.debug_default is not None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += "{} = {};\n".format(tstr, self.debug_default)
            ostr += "#else\n"
            ostr += "{} = {};\n".format(tstr, self.default)
            ostr += "#endif\n"
        else:
            if self.dtype == "string" and self.default == '""':
                ostr += "{};\n".format(tstr)
            elif self.dtype == "dim_array":
                ostr += "{} = {{{}}};\n".format(tstr, self.default)
            else:
                ostr += "{} = {};\n".format(tstr, self.default)

        if self.ifdef is not None:
            ostr += "#endif\n"

        return ostr

    def get_query_string(self, language):
        # this is the line that queries the ParmParse object to get
        # the value of the runtime parameter from the inputs file.
        # This goes into pelec_queries.H included into PeleC.cpp

        ostr = ""
        if self.ifdef is not None:
            ostr += "#ifdef {}\n".format(self.ifdef)

        if language == "C++":
            if self.dtype == "dim_array":
                ostr += "{\n  amrex::Vector<amrex::Real> tmp(AMREX_SPACEDIM, 0.0);\n"
                ostr += '  pp.queryarr("{}", tmp, 0, AMREX_SPACEDIM);\n'.format(
                    self.name
                )
                ostr += "  for (int i = 0; i < tmp.size(); i++) {{\n    {}[i] = tmp[i];\n  }}\n".format(
                    self.cpp_var_name
                )
                ostr += "}\n"
            else:
                ostr += 'pp.query("{}", {});\n'.format(self.name, self.cpp_var_name)
        else:
            sys.exit("invalid language choice in get_query_string")

        if self.ifdef is not None:
            ostr += "#endif\n".format(self.ifdef)

        return ostr

    def get_decl_string(self):
        # this is the line that goes into pelec_params.H included
        # into PeleC.H

        static = ""
        if self.static:
            static = "static"

        if self.dtype == "int":
            tstr = "{} int {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "{} amrex::Real {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "{} std::string {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "{} bool {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "dim_array":
            tstr = "{} amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> {};\n".format(
                static, self.cpp_var_name
            )
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if self.ifdef is not None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        ostr += tstr

        if self.ifdef is not None:
            ostr += "#endif\n"

        return ostr


def parse_params(infile):

    params = []

    namespace = None
    cpp_class = None
    static = None

    try:
        f = open(infile)
    except IOError:
        sys.exit("error opening the input file")

    for line in f:
        if line[0] == "#":
            continue

        if line.strip() == "":
            continue

        if line[0] == "@":
            # this is a command
            cmd, value = line.split(":")
            if cmd == "@namespace":
                fields = value.split()
                namespace = fields[0]
                cpp_class = fields[1]

                try:
                    static = fields[2]
                except IndexError:
                    static = ""

                # do we have the static keyword?
                if "static" in static:
                    static = 1
                else:
                    static = 0

            else:
                sys.exit("invalid command")

            continue

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r"[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)", line)

        name = fields[0]
        if name[0] == "(":
            name, cpp_var_name = re.findall(r"\w+", name)
        else:
            cpp_var_name = name

        dtype = fields[1]

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try:
            ifdef = fields[4]
        except IndexError:
            ifdef = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(
            Param(
                name,
                dtype,
                default,
                cpp_var_name=cpp_var_name,
                namespace=namespace,
                cpp_class=cpp_class,
                static=static,
                debug_default=debug_default,
                ifdef=ifdef,
            )
        )

    # output

    # find all the namespaces
    namespaces = list(set([q.namespace for q in params]))

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]

        # write name_defaults.H
        try:
            cd = open("{}/{}_defaults.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_defaults.H for writing".format(nm))

        cd.write(CWARNING)

        for p in params_nm:
            cd.write(p.get_default_string())

        cd.close()

        # write name_params.H
        try:
            cp = open("{}/{}_params.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_params.H for writing".format(nm))

        cp.write(CWARNING)

        for p in params_nm:
            cp.write(p.get_decl_string())

        cp.close()

        # write pelec_queries.H
        try:
            cq = open("{}/{}_queries.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_queries.H for writing".format(nm))

        cq.write(CWARNING)

        for p in params_nm:
            cq.write(p.get_query_string("C++"))

        cq.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        nargs=1,
        help="input file containing the list of parameters we will define",
    )
    args = parser.parse_args()

    parse_params(args.input_file[0])
