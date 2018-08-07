""" Freeing up space for new observations by deleting old files. """


import ctypes
import os
import platform
import sys
import shutil
import datetime
import math
import ephem

from CaptureDuration import captureDuration



def availableSpace(dirname):
    """
    Returns the number of free bytes on the drive that p is on.

    Source: https://atlee.ca/blog/posts/blog20080223getting-free-diskspace-in-python.html
    """

    if platform.system() == 'Windows':

        free_bytes = ctypes.c_ulonglong(0)

        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(dirname), None, None, \
            ctypes.pointer(free_bytes))

        return free_bytes.value

    else:
        st = os.statvfs(dirname)

        return st.f_bavail*st.f_frsize




def getNightDirs(dir_path):
    """ Returns a sorted list of directories in the given directory which conform to the captured directories
        names. 

    Arguments:
        dir_path: [str] Path to the data directory.
        stationID: [str] Name of the station. The directory will have to contain this string to be taken
            as the night directory.

    Return:
        dir_list: [list] A list of night directories in the data directory.

    """

    # Get a list of directories in the given directory
    dir_list = [dir_name for dir_name in os.listdir(dir_path)]
    dir_list = sorted(dir_list)

    return dir_list



def deleteNightFolders(dir_path, delete_all=False):
    """ Deletes captured data directories to free up disk space. Either only one directory will be deleted
        (the oldest one), or all directories will be deleted (if delete_all = True).

    Arguments:
        dir_path: [str] Path to the data directory.

    Keyword arguments:
        delete_all: [bool] If True, all data folders will be deleted. False by default.

    Return:
        dir_list: [list] A list of remaining night directories in the data directory.

    """

    # Get the list of night directories
    dir_list = getNightDirs(dir_path)

    # Delete the saving directories if the macro file, CaptureData or ArchivingData
    for dir_name in dir_list:
        bytes_before = availableSpace(data_dir)
        shutil.rmtree(os.path.join(dir_path,dir_name))
        byte_change = bytes_before - availableSpace(data_dir)
        print("Removing {s} to free {:} bytes, leaving {:} free bytes.".format(os.path.join(dir_path,dir_name), byte_change, availableSpace(data_dir)))
        # If only one (first) file should be deleted, break the loop
        if not delete_all:
            break

    # Return the list of remaining night directories
    return getNightDirs(dir_path)



def deleteOldObservations(data_dir, captured_dir, archived_dir, duration):
    """ Deletes old observation directories to free up space for new ones.

    Arguments:
        data_dir: [str] Path to the RMS data directory which contains the Captured and Archived diretories
        captured_dir: [str] Captured directory name.
        archived_dir: [str] Archived directory name.
        duration: [float] Duration of next video capturing in seconds. If None (by default), duration will
            be calculated for the next night.

    Return:
        [bool]: True if there's enough space for the next night's data, False if not.

    """

    captured_dir = os.path.join(data_dir, captured_dir)
    archived_dir = os.path.join(data_dir, archived_dir)
    # Bytes per file
    # Estimate the size of the file by summing all components of the rdm struct
    RDM_SIZE = (32 + 32 + 6*8 + 1*8 + 64 + 64 + 64 + 64*8 + 32 + 64 + 32 + 32 + 32 + 32 + 32*(1024**2) + 32*(1024**2) + 32*(1024**2))/8;
    
    # Calculate the number of files that will be saved tonight
    # 7.2 files per hour 3600(s/hr)/(500s/file) = 7.2 files/hr 
    loops = int(math.ceil(7.2*duration))
    
    # Calculate the approx. size for the night 
    next_night_bytes = int(loops*RDM_SIZE)

    # Always leave at least 1 GB free for archive
    next_night_bytes += 5*(1024**3)
    
    print("Space required for the next night is {:} bytes\nand the space available is {:} bytes.".format(next_night_bytes, availableSpace(data_dir)))

    # If there's enough free space, don't do anything
    if availableSpace(data_dir) > next_night_bytes:
        print("There is adequate space available.")
        return True
    
    # Get the current directories
    captured_dirs_remaining = getNightDirs(captured_dir)

    # Intermittently delete captured and archived directories until there's enough free space
    while True:

        # If there is still another directory to delete 
        if not (len(captured_dirs_remaining) == 0):
            
            # Delete one captured file or an empty captured directory
            captured_dirs_remaining = deleteNightFolders(captured_dir)
            archived_dirs_remaining = deleteNightFolders(archived_dir)
            
            # Break the there's enough space
            if availableSpace(data_dir) > next_night_bytes:
                break

            if ((len(captured_dirs_remaining) == 0) or (len(archived_dirs_remaining) == 0)):
                print("Program terminating, insufficient memory to run program. Memory must be cleared manually.")
                sys.exit()
            
                
    return True

if __name__ == "__main__":
    path = "/home/pi/RadiometerData"
    capt = "CapturedData"
    arch = "ArchivedData"
    deleteOldObservations(path,capt, arch, 4)
