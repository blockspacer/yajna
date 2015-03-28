#!/usr/bin/python

import itertools
import os
import fnmatch
import re
import sys
import argparse

import ninja_syntax
import gcc
import msvc

# --- util functions

def flags(*iterables):
    return ' '.join(itertools.chain(*iterables))

def get_files(root, pattern):
    pattern = fnmatch.translate(pattern)
    for dir, dirs, files in os.walk(root):
        for f in files:
            if re.match(pattern, f):
                yield os.path.join(dir, f)

def object_file(fn):
    return os.path.join('obj', re.sub(r'\.c\+\+$', '.o', fn))

# --- arguments

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true', help='compile with debug information')
parser.add_argument('--cxx', default=None, metavar='executable', help='compiler name to use')
parser.add_argument('--msvc', action='store_true', help='use the MSVC++ toolchain')
parser.add_argument('--boost-dir', default=None, metavar='path', help='path of boost folder (i.e. the folder with include/ and lib/ subfolders)')
parser.add_argument('--no-lto', action='store_true', help='do not perform link-time optimisation')
args = parser.parse_args()

tools = msvc.Toolchain() if args.msvc else gcc.Toolchain()
compiler = args.cxx if args.cxx else tools.compiler()
linker = args.cxx if args.cxx else tools.linker()

# --- variables

dependencies = []
include_flags = flags([tools.include('include')],
                      (tools.dependency_include(os.path.join('deps', d, 'include')) for d in dependencies))
if(args.boost_dir):
    include_flags += ' ' + tools.dependency_include(args.boost_dir)
cxx_flags = flags(tools.cxx_flags(),
                  tools.debug_flags() if args.debug else tools.optimisation_flags(),
                  [] if args.no_lto or args.debug else tools.linker_lto_flags())
warning_flags = flags(tools.max_warnings())
define_flags = ''
lib_flags = ''
ld_flags = flags(tools.link_flags(),
                 [] if args.no_lto or args.debug else tools.linker_lto_flags())

stringize_tool = 'tools/stringize.py' 
single_header_tool = 'tools/single_header.py'

# --- preamble

ninja = ninja_syntax.Writer(open('build.ninja', 'w'))

ninja.variable('ninja_required_version', '1.3')
ninja.variable('builddir', 'obj' + os.sep)
ninja.variable('msvc_deps_prefix', 'Note: including file:')
 
# --- rules

ninja.rule('bootstrap',
        command = ' '.join(['python'] + sys.argv),
        generator = True,
        description = 'BOOTSTRAP')

ninja.rule('cxx',
        command = ' '.join([compiler, flags(tools.dependencies_output('$out.d')), cxx_flags, warning_flags, include_flags, define_flags, '$extraflags', '$in', flags(tools.compiler_output('$out'))]),
        deps = tools.ninja_deps_style(),
        depfile = '$out.d',
        description = 'C++ $in')

ninja.rule('link',
        command = ' '.join([linker, ld_flags, '$in', flags(tools.linker_output('$out')), lib_flags]),
        description = 'LINK $out')

# --- build edges

ninja.build('build.ninja', 'bootstrap',
        implicit = sys.argv[0])

hdr_files = list(get_files('include', '*.h++'))
src_files = list(get_files('src', '*.c++'))
obj_files = [object_file(fn) for fn in src_files]
for fn in src_files:
    ninja.build(object_file(fn), 'cxx',
            inputs = fn)

program = os.path.join('bin', 'yajna') + tools.executable_extension()
ninja.build(program, 'link',
        inputs = obj_files)
ninja.build('yajna', 'phony',
        inputs = program)

test_src_files = list(get_files('test', '*.c++'))
test_obj_files = [object_file(fn) for fn in test_src_files]
for fn in test_src_files:
    ninja.build(object_file(fn), 'cxx',
            inputs = fn)

test_runner = os.path.join('bin', 'test') + tools.executable_extension()
ninja.build(test_runner, 'link',
        inputs = test_obj_files)
ninja.build('test', 'phony',
        inputs = test_runner)

ninja.default('yajna')

