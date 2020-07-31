#!/usr/bin/env python

import os, sys, shutil, stat
import argparse, subprocess

# Credit for this function to George King (https://github.com/gwk)
# Lifted from https://github.com/gwk/gloss/blob/master/python/gloss/otool.py

# Edit history:
# 2020-07-31: Started keeping edit history
# 2020-07-31: Added error reporting for

def otool(s):
    o = subprocess.Popen(['otool', '-L', s], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error1 = ''.join(i for i in o.stderr.readline())
    if not "" == error1:
        print("Error:" + error1)
    for l in o.stdout:
        if l[0] == '\t':
            yield l.split(' ', 1)[0][1:]


def getLibList(lib, cleaner=None):
    done=set()
    left=set([lib])
    while not len(left) == 0:
        tLib = left.pop()
        done.add(tLib)
        print "Inspecting " + tLib
        s = set(otool(tLib))
        if cleaner is not None:
            s = cleaner(s)
        print "     Got: " + str(s)
        s -= done
        left.update(s)
    return done

def consolidateLibs(libs, tdir=["/usr/local/lib"], sdir=".", link_path="@loader_path"):
    moved = []
    for lib in libs:
        folder, fname = os.path.split(lib)
        for comp in tdir:
            if folder.startswith(comp):
                shutil.copy(os.path.join(folder, fname), os.path.join(sdir, fname))
                moved += [[os.path.join(folder, fname), os.path.join(link_path, fname)]]
                os.chmod(os.path.join(sdir, fname), stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
                next
    print "Moved: " + str(moved)
    return moved

def updateLibs(libs, moved):
    for source, target in moved:
        for tlib in libs:
            lib = os.path.basename(tlib)
            print "Updating " + lib + " to reflect move from " + source + " to " + target
            pipe=subprocess.Popen(['install_name_tool', '-change', source, target, lib],
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True)
            error1 = ''.join(i for i in pipe.stderr.readline())
            if not "" == error1:
                print("Error:" + error1)


def updateIDs(libs):
    for source, tlib in libs:
        lib = os.path.basename(tlib)
        print "Updating " + source + " with name " + lib
        pipe = subprocess.Popen(['install_name_tool', '-is', lib, lib],
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True)
        error1 = ''.join(i for i in pipe.stderr.readline())
        if not "" == error1:
            print("Error:" + error1)

def make_cleaner(bads):
    def cleaner(names):
        # print "Called cleaner."
        outNames = set()
        for name in names:
            for bad in bads:
                # print( " Does " + name + " start with "  + bad + "? " + str(name.startswith(bad)))
                if name.startswith(bad):
                    outNames.add(name)
        return outNames
    return cleaner

if __name__ == "__main__":
    print "Welcome to libtastic"
    parser = argparse.ArgumentParser(description="Libtastic traces OS X libraries and moves and adjusts them.", usage="libtastic.py -id OpenMx.so <library_root_dir>",epilog="")
    parser.add_argument('-id', '--updateIDs', '--updateids', action="store_true")
    parser.add_argument('lib', help='Initial Library')
    parser.add_argument('locs', nargs='*', default="/opt/local/lib/")

    args = parser.parse_args()

    print("lib=" + args.lib)
    print("locs=" + str(args.locs))
    libList = getLibList(args.lib, make_cleaner(args.locs))
    # moved = consolidateLibs(libList, args.locs)
    updateLibs(libList, zip(libList, libList))# moved)
    updateIDs(moved)
    print "Thank you for being libtastic"
