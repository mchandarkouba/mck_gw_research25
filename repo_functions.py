###############################################################################
### IMPORTS
###############################################################################

import os
import datetime
from git import Repo # install gitpython to env

###############################################################################
### GLOBALS
###############################################################################

PROJ_DIR = "/Users/mck/Desktop/palmese_research/"
REPO_DIR = "/Users/mck/Desktop/palmese_research/research_repo/"
DATE = datetime.datetime.now().strftime( "%m-%d-%Y" )

###############################################################################
### HELPERS
###############################################################################

def recursive_all_contents(ldir, 
                           handling_type="all", 
                           exclude_dirs=None,
                           ):
    result = []
    for s in os.listdir(ldir):
        path = ldir+s
        
        add_to_result = handling_type=="files" and os.path.isfile(path)
        add_to_result |= handling_type=="dirs" and os.path.isdir(path)
        add_to_result |= handling_type=="all"
        
        if add_to_result: result += [path]
        
        if os.path.isdir(path):
            if isinstance(exclude_dirs, str) and os.path.samefile(path,exclude_dirs): continue
            elif isinstance(exclude_dirs, list) and any( os.path.samefile(path, x) for x in exclude_dirs): continue
            result += recursive_all_contents(path+'/')
            
    return result

###############################################################################

def namesearch_path(name,
                    directory,
                    handling_type="all",
                    exclude_dirs=None,
                    ):
    if name==None: raise NameError
    
    dirlist = recursive_all_contents(directory,
                                     handling_type=handling_type,
                                     exclude_dirs=exclude_dirs,
                                     )
    results = list()
    
    for file in dirlist: 
        if name in file: results.append(file)

    if len(results)==1: return results[0]
    
    elif len(results)==0:
        print(f"No files were found for name ({name}) search")
        raise NameError
    
    else:
        return choose_path(results)

###############################################################################

def choose_path(dirlist):
    master = {n+1:p for n,p in enumerate(dirlist)}
    
    in_bounds = lambda n: n.isnumeric() and int(n)>=1 and int(n)<=len(dirlist)
    print("Choose preferred file (by index number)")
    for n in master: print(f"\t{n} : {master[n]}")
    
    c = 0
    while True:
        c=input('index: ')
        if in_bounds(c): 
            c=int(c)
            break
    
    print(f"\nUsing {master[c]}")
    return master[c]

###############################################################################

def update_changelog(path_in=None, 
                     path_out=None,
                     add_note='',
                     mode:str=None,
                     ):
    generic = ""    
    hbar = '\n' + '-'*50 + '\n'

    if mode=="replace":
        generic = f"Copied:\n\t{path_in}\n\t--> {path_out}\n"
        
    elif mode=="add":
        generic = f"Added:\n\t{path_in}"
    
    elif mode=="remove":
        generic = f"Removed:\n\t{path_in}"
    
    dir_out = os.path.dirname(path_out)
    changelog = dir_out + "/changelog.txt"
    if not os.path.exists(changelog): os.makedirs(changelog)
    
    with open(changelog, 'r') as file: lines = file.readlines()
    
    new_update = True
    
    if len(lines)>=3:
        rlines = lines[::-1]
        
        for i, in range(len(rlines)-2):
            l0, l1, l2 = rlines[i], rlines[i+1], rlines[i+2]
            if l0==l2==hbar and l1==DATE: new_update = False
    
    update = hbar+DATE+hbar+generic if new_update else generic
    update = update + '\n\n' + add_note + '\n'
    with open(changelog, 'a') as file: file.write(update)

###############################################################################

def verify_process():
    while True:
        verification = input("Proceed? (y/n) ")
        
        if verification =='n':
            print("Cancelling process")
            raise Exception
            
        if verification == 'y': return

###############################################################################
### FUNCTIONS
###############################################################################

def copy_to_repo(name_in=None, name_out=None, name=None, 
                 path_in=None, path_out=None,
                 log_note='',
                 verify=True,
                 handling_type="all",
                 ):
    
    if name_in==None and name_out==None: name_in = name_out = name
    
    if path_in==None: path_in = namesearch_path(name_in, 
                                                PROJ_DIR, 
                                                handling_type=handling_type, 
                                                exclude_dirs=REPO_DIR,
                                                )
    if path_out==None: path_out = namesearch_path(name_out,
                                                  REPO_DIR,
                                                  handling_type=handling_type,
                                                  )
    
    print(f"Attempt to copy:\n\t{path_in}\n\t--> {path_out}\n")
    if verify: verify_process()
    
    flag = "-r " if os.path.isdir(path_in) else ''
        
    try:
        os.system(f"cp {flag}{path_in} {path_out}")
        
    except Exception as e: 
        print(f"Error during copy: {e}")
    
    update_changelog(path_in=path_in,
                     path_out=path_out,
                     add_note=log_note, 
                     mode="replace",
                     )
        
###############################################################################

def add_to_repo(path=None,
                 log_note=None,
                 verify=True,
                 ):
    print(f"Attempt to add:\n\t{path}")
    if verify: verify_process()
    
    try:
        os.system(f"touch {path}")
        
    except Exception as e: 
        print(f"Error during creation: {e}")
        
    update_changelog(path_in=path,
                     add_note=log_note, 
                     mode="add",
                     )

###############################################################################

def remove_from_repo(name=None, path=None,
                     log_note="",
                     verify=True,
                     handling_type="all",
                     ):
    if path==None: path = namesearch_path(name, 
                                          REPO_DIR, 
                                          handling_type=handling_type, 
                                          )
    print(f"Attempt to remove:\n\t{path}")
    if verify: verify_process()
    
    command="rmdir" if os.path.isdir(path) else "rm"
    
    try:
        os.system(f"{command} {path}")
        
    except Exception as e: 
        print(f"Error during deletion: {e}")
        
    update_changelog(path_in=path,
                     add_note=log_note, 
                     mode="remove",
                     )

###############################################################################

def push_github(message=None):
    generic = f"Update {DATE}"
    message = message or generic
    
    try:
        repo = Repo(REPO_DIR)
        repo.git.add(A=True)
        repo.index.commit(message)
        origin = repo.remote(name="origin")
        origin.push(refspec="main:main")
        
    except Exception as e:
        print(f"Error pushing to origin: {e}")
    
    return

###############################################################################

def tests():
    global REPO_DIR
    REPO_DIR = "/Users/mck/Desktop/palmese_research/test_repo/"
    
    copy_to_repo(name_in="observation_populations", path_out=REPO_DIR+"copy_dir", handling_type="all")
