#!/usr/bin/env python

import os, sys, shutil, stat
import argparse, subprocess, re
import errno
import traceback
import json

# Credit for this function to George King (https://github.com/gwk)
# Lifted from https://github.com/gwk/gloss/blob/master/python/gloss/otool.py

# Edit history:
# 2020-07-31: Started keeping edit history
# 2020-07-31: Added error reporting for otool and install_name_tool calls

def getRpaths(s, root=None, rpaths=[]):
    libPaths = rpaths[:]
    localPath = os.path.dirname(s)
    if localPath:
        libPaths.append(localPath)
    if root is not None:
        libPaths.append(os.path.dirname(root))
    libPaths.append('.')
    libPaths.append(os.path.join(localPath, '..'))
    o = subprocess.Popen(['otool', '-l', s], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error1 = ''.join(i for i in o.stderr.readline())
    if not "" == error1:
        print("Error:" + error1)
        print("Otool called on " + s + " from:")
        print(traceback.print_stack())
    foundRpath = False
    
    # Matches the mac RPATH spec.  Specifically:
    #      cmd LC_RPATH
    #  cmdsize ##
    #     path [THE THING WE CARE ABOUT] (offset ##)\n
    # And does it in a way that ignores spacing and won't bleed between rpath specs
    pathmatch = re.compile("cmd LC_RPATH\n[^\n]*\n[\s]*path ([\s\S]*?) \(offset [\w]+?\)\n")
    rpath_hit = re.compile("@loader_path|@executable_path|@rpath")
    rpath_miss = re.compile("((?!@loader_path|@executable_path|@rpath).*)")
    libPaths += [x for x in pathmatch.findall(o.stdout.read())]
    # print(libPaths)
    for aPath in filter(rpath_hit.match, libPaths):
        for bPath in filter(rpath_miss.match, libPaths):
            newPath = rpath_hit.sub(bPath, aPath, 1)
            if newPath != aPath:
                libPaths.append(newPath)
                
    # print("Rpaths for " + s + " are " + str(libPaths))
    return libPaths

def otool(s):
    o = subprocess.Popen(['otool', '-L', s], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error1 = ''.join(i for i in o.stderr.readline())
    if not "" == error1:
        print("Error:" + error1)
        print("Otool called on " + s + " from:")
        print(traceback.print_stack())
    output = []
    for l in o.stdout:
        if l[0] == '\t':
            output.append(l.split(' ', 1)[0][1:])
    return output[1:]

def getLibList(lib, keepOnly=None, *args, **kwargs):
    ''' otool through the tree to get name/target pairs.
        If requested, the tree can be pruned by keeping only
        libraries that fall in a certain set of folders.
        
        Output is a list of dictionaries with three
        elements: 
            short (the name as displayed in the library)
            full  (the resolved path to a real file)
            root  (the resolved path of the library from which the file is linked)
    '''

    resolved = []
    completed=[]
    left=[lib]
    sources=['.']
    rpath_regex = re.compile("(@loader_path|@executable_path|@rpath)")

    while not len(left) == 0:
        tLib = left.pop()
        print "Inspecting " + str(tLib)
        if tLib in completed:  # Did this one already
            next
        needed = otool(tLib) # Does the actual work
        # print("Linked libraries for " + tLib + " are " + str(needed))
        possibles = getRpaths(tLib, root=sources.pop(), *args, **kwargs)
        # print("Possible resolutions for " + tLib + " are " + str(possibles))
        for aPath in needed:
            if os.path.isfile(aPath):
                print("   Found " + aPath + " in place.")
                found=True
                selected=aPath
            elif rpath_regex.match(aPath):
                selected = None
                found = False
                for tgt in possibles:
                    # print("   Replacing " + rpath_regex.match(aPath).group(0) + " with " + tgt + " in " + aPath + " to get " + rpath_regex.sub(tgt, aPath, 1))
                    newPath = rpath_regex.sub(tgt, aPath, 1)
                    # print("  Trying: " + newPath)
                    if os.path.isfile(newPath):
                        print("  Found: " + newPath + " by replacing " + rpath_regex.match(aPath).group(0) + " with " + tgt)
                        found = True
                        selected = newPath
                        break
            else:
                print("Library " + str(aPath) + " does not resolve.")
                raise Exception("All resolutions of " + str(aPath) + " for file " + tLib + " failed.")
            if not found:
                print("I tried rpaths " + str(possibles) + " for file " + tLib + ".")
                raise Exception("All resolutions of " + str(aPath) + " for file " + tLib + " failed.")
            if selected is not None:   # Handle cases where only keeping specific ones.
                # print("Resolving " + aPath + " as " + selected)
                resolved.append({'short':aPath, 'full':selected, 'root':tLib})
                if keepOnly:
                    for keeper in keepOnly:
                        if selected.startswith(keeper):
                            left.insert(0,selected)
                            sources.insert(0,tLib)
                            print("    Also searching " + selected)
                else:
                    left.insert(0,selected)
                    sources.insert(0,tLib)
            else:
                print("Nothing selected for " + aPath)
                
    return resolved
    
def showLibList(lib, cleaner=None):
    fullPath = os.path.abspath(lib)
    print("Library list for library at " + fullPath)
    print(set(otool(fullPath)))

def touch(fname, times=None):
        with open(fname, 'a'):
            os.utime(fname, times)
    
def consolidateLibs(libs, target={'source':"/usr/local/lib", 'target':".",  'link_path':"@loader_path"}):
    ''' Moves libraries.
    
        libs is the output of getLibList.
    
        target provides a dict:
            source: any library in this folder will be moved
            target: to here (usually '.' or the path of the loading library)
            link_path: and links updated with this path (usually @loader_path)
    
        returns libs, updated with moved and link paths.
    '''
    newLibs = []
    if not libs:   # Nothing to do.
        print("No libraries to consolidate.  Skipping.")
        return libs
    rootmap = {}
    for lib in libs:
        if lib["root"] in rootmap.keys():  # adjust if the root has been moved.
            lib["root"] = rootmap[lib["root"]]
        folder, fname = os.path.split(lib['full'])
        for link in target:
            
            if folder.startswith(link["source"]):
                print("Copying:" + str((os.path.join(folder, fname), os.path.join(link["target"], fname))))
                shutil.copy(os.path.join(folder, fname), os.path.join(link["target"], fname))
                lib["moved"] = os.path.join(link["target"], fname)
                lib["link_path"] = os.path.join(link['link_path'], fname)
                os.chmod(os.path.join(link["target"], fname), stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)
                print "Moved: " + lib['full'] + " to " + lib['moved'] + " (called " + lib['short'] + " in library " + lib['root'] + ") to link as " + lib['link_path']
                newLibs.append(lib)
                rootmap[lib['full']] = lib['moved']
    print("Consolidated: " + str(newLibs))
    return newLibs

def updateLibs(libs):
    ''' Uses install_name_tool to update the internal links and names
        of the library to reflect the moves.
    
        input is the list of libraries, old and new locations, and link paths
        coming from consolidateLibs.
    '''
    if not libs:
        print("No libraries to update.  Skipping.")
        return
    for tlib in libs:
        print("Updating "+ str(tlib))
        lib = os.path.basename(tlib['short'])
        print "Updating " + tlib['short'] + " to " + tlib['link_path'] + " in " + tlib["root"] + " to reflect move from " + tlib['full'] + " to " + tlib['moved']
        pipe=subprocess.Popen(['install_name_tool', '-change', tlib['short'], tlib['link_path'], tlib['root']],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
        error1 = ''.join(i for i in pipe.stderr.readline())
        if not "" == error1:
            print("Error:" + error1)

def updateIDs(libs):
    if not libs:
        print("No IDs to update. Skipping.")
        return
    for tlib in libs:
        lib = os.path.basename(tlib['short'])
        src = tlib['moved']
        print "Updating " + tlib['moved'] + " with name " + lib
        pipe = subprocess.Popen(['install_name_tool', '-id', lib, tlib['moved']],
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True)
        error1 = ''.join(i for i in pipe.stderr.readline())
        if not "" == error1:
            print("Error:" + error1)


if __name__ == "__main__":
  try:
    print "Welcome to libtastic"
    parser = argparse.ArgumentParser(description="Libtastic traces OS X libraries and moves and adjusts them.", usage="libtastic.py -id OpenMx.so <library_root_dir>",epilog="")
    parser.add_argument('-id', '--updateIDs', '--updateids', action="store_true", help="update library links to new locations (omit for dry run)")
    parser.add_argument('-rp', '--rpath', help="Search list of folders to follow for @rpath links", action ="append", default=[str(os.environ["HOME"])+"/lib:/usr/local/lib:/lib:/usr/lib"])
    parser.add_argument('lib', help='Initial Library')
    parser.add_argument('locs', nargs='*', default="/opt/local/lib")

    args = parser.parse_args()
    
    args.locs = [x.rstrip("/") for x in args.locs] # kill trailing /
    print("lib=" + str(args.lib))
    print("locs=" + str(args.locs))
    locList = [{'source':loc, 'target':'.', 'link_path':'@loader_path'} for loc in args.locs]
    rpaths = [x for y in args.rpath for x in y.split(":") ]
    rpaths += args.locs
    print("rpaths=" + ", ".join(rpaths))
    showLibList(args.lib)
    libList = getLibList(args.lib, keepOnly=args.locs, rpaths=rpaths)
    # print("Processed.  Got: \n" + json.dumps(libList) )
    # print("Consolidating for " + str(locList))
    moved = consolidateLibs(libList, locList)
    updateLibs(moved)
    # print(moved)
    if(args.updateIDs):
        updateIDs(moved)
    showLibList(args.lib)
    touch("EditedByLibtastic.txt")
    print "Thank you for being libtastic"
  except Exception:
    print(traceback.format_exc())
