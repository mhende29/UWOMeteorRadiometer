
import os
import shutil
import tarfile
import time

def manArchive(source_dir, dest_dir):
    
    file_list = os.listdir(source_dir)
    
    files_zipped = 0
    plots_moved = 0
    
    rdm_list = [rdm_file for rdm_file in file_list if rdm_file.endswith(".rdm")]
    rdm_list = sorted(rdm_list)
    plot_list = [plot for plot in file_list if plot.endswith(".png")]
    
    for rdm_file in rdm_list:
        # Compress the archive 
        archive_file_name = os.path.join(source_dir, rdm_file + ".tar.bz2")
        
        with tarfile.open(archive_file_name, "w:bz2") as tar:
            tar.add(os.path.join(source_dir,rdm_file))
        
        shutil.move(archive_file_name,os.path.join(dest_dir, os.path.split(archive_file_name)[1]))
        
        files_zipped += 1
        
        if(rdm_file is not rdm_list[-1]):
            print("Zipped {:.2%} of files.".format(files_zipped/len(rdm_list)), end = '\r')
        else:
            print("Zipped {:.2%} of files.".format(files_zipped/len(rdm_list)), end = '\n')

    for plot in plot_list:
        
        shutil.move(os.path.join(source_dir, plot), os.path.join(dest_dir, plot))
        
        plots_moved += 1        
        
        if(plot is not plot_list[-1]):
            print("Moved {:.2%} of plots.".format(plots_moved/len(plot_list)), end = '\r')
        else:
            print("Moved {:.2%} of plots.".format(plots_moved/len(plot_list)), end = '\n')





if __name__ == "__main__":
    source1 = "/home/pi/RadiometerData/CapturedData/CA0001_A_20180707-014143"
    source2 = "/home/pi/RadiometerData/CapturedData/CA0001_A_20180708-014116"
    source3 = "/home/pi/RadiometerData/CapturedData/CA0001_A_20180709-014047"
    
    dest1 = "/home/pi/RadiometerData/ArchivedData/CA0001_A_20180707-014143"
    dest2 = "/home/pi/RadiometerData/ArchivedData/CA0001_A_20180708-014116"
    dest3 = "/home/pi/RadiometerData/ArchivedData/CA0001_A_20180709-014047"
    
    #manArchive(source1,dest1)
    #manArchive(source2,dest2)
    manArchive(source3,dest3)
