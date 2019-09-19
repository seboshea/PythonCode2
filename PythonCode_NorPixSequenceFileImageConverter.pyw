import os
import datetime
import tkinter as tk
import tkinter.filedialog
import numpy as np
#import logging
import gettext
import multiprocessing
import functools
import PIL.Image
import PIL.PngImagePlugin


#NorPixSequenceFileImageConverter(sequenceFilePath='C:/Users/Admin TEMP/Documents/HaloHolo/Test/', imageDir='C:/Users/Admin TEMP/Documents/HaloHolo/Test/Output/', imagePrefix='TEST_', overWrite=False, cropcoords=[])
def NorPixSequenceFileImageConverter(sequenceFilePath='', imageDir='', imagePrefix='C075_', overWrite=False, cropcoords=[]):
    """ A small program to convert .seq files, taken by HALOHOLO, to .png holorams"""
    
    if len(cropcoords) not in [0,4]:
        print('Cropcoords have not the right number (4) of values. Falling back to default ([]).')
        cropcoords = []
    if len(sequenceFilePath) == 0 or not os.path.isdir(sequenceFilePath):
        seqDir = tk.filedialog.askdirectory()
        if len(seqDir) == 0:
            print('User failed to choose a valid folder for *.seq files\n Operation cancelled by user')
            return
    else:
        seqDir = sequenceFilePath
    sequenceFilePaths = [file for file in os.listdir(os.path.normpath(os.path.join(seqDir))) if file.endswith('.seq')]
    if len(imageDir) == 0 or not os.path.isdir(imageDir):
        imageDir = tk.filedialog.askdirectory()
        if len(imageDir) == 0:
            print('User failed to choose a valid folder for holograms production\n Operation cancelled by user')
            return

    isfinished = False

    while not isfinished:
        isfinished = processSeqFiles(seqDir,sequenceFilePaths, imageDir, imagePrefix, overWrite, cropcoords)


def processSeqFiles(seqDir, sequenceFilePaths, imageDir, imagePrefix, overWrite, cropcoords):
    doneSeqFiles = loadBlacklist(imageDir)
    sequenceFilePaths = [x for x in sequenceFilePaths if x not in doneSeqFiles]

    # Fill out to full search paths
    sequenceFilePaths = [os.path.join(seqDir,x) for x in sequenceFilePaths]

    nSeqFiles = len(sequenceFilePaths)
    if nSeqFiles == 0:
        return True
#   If multiprocessing is not wanted use this and not the following code        
#    for file in sequenceFilePaths:
#        processOneSequenceFile(file, imageDir, imagePrefix, overWrite, cropcoords)

    with multiprocessing.Pool() as pool:
        pool.map(functools.partial(processOneSequenceFile,imageDir=imageDir, imagePrefix=imagePrefix, overWrite=overWrite, cropcoords=cropcoords),sequenceFilePaths)

    return False
        

def processOneSequenceFile(thisSequenceFilePath, imageDir, imagePrefix, overWrite, cropcoords):
    fileSize = os.path.getsize(thisSequenceFilePath)
    
    with open(thisSequenceFilePath,'rb') as fid:
        MagicNumber = np.fromfile(fid,dtype='<u4',count=1,sep='')[0]
        if MagicNumber!=65261:
            print("{} has the wrong magic number for a sequence file".format(thisSequenceFilePath))
            return 1
#        fid.seek(32, 0)
        fileHeaderLength = 8192

        fid.seek(548,0)
        Magic = np.fromfile(fid,dtype='<u4',count=6,sep='')
        imageSizeNx = Magic[0]
        imageSizeNy = Magic[1]
        imageBitDepth = Magic[2]
        imageBitDepthReal = Magic[3]
        imageLengthBytes = Magic[4]
        imageFormat = Magic[5]

        fid.seek(572, 0)
        # Get the number of images that are "supposed" to be in this file
        numImages = np.fromfile(fid,dtype='<u4',count=1,sep='')[0]

        if (imageBitDepth != 8)|(imageBitDepthReal != 8)|(imageFormat != 100):
            print('This function only supports Monochrome 8-bit images')

        fid.seek(580,0)
        # Get the spacing between images as they are separated by the image size
        # + the image footer up to the next sector boundary
        imageSpacing = np.fromfile(fid,dtype='<u4',count=1,sep='')[0]

        # Suggested frame rate
#        fid.seek(584,0)
#        frameRate = fid.read(1).astype('double')

        # Suggested frame rate ? extended Header
#        extendedHeader = fid.read(1,0).astype('int32')

        time_t_roottime = datetime.datetime(1970,1,1,0,0,0) + datetime.timedelta(days=366)

        shouldBeNumImages = int((fileSize - fileHeaderLength)/(imageSpacing))
        if numImages != shouldBeNumImages:
            print('File may be corrupted, the number of reported images does not match the file size.')

        numImages = shouldBeNumImages

#        logger = logging.getLogger(__name__)

        for cnt2 in range(numImages):
            offset = fileHeaderLength + cnt2*imageSpacing
            fid.seek(offset+imageLengthBytes,0)
            rawImTime1 = np.fromfile(fid,dtype='<u4',count=1,sep='')[0]
            rawImTime2 = np.fromfile(fid,dtype='<u2',count=1,sep='')[0]
            rawImTime3 = np.fromfile(fid,dtype='<u2',count=1,sep='')[0]
            rawImTime = rawImTime1 + rawImTime2/1000 + rawImTime3/1e6
            imTime = rawImTime #+ time_t_roottime.toordinal()
            imFilename = "{}{}.png".format(imagePrefix,datetime.datetime.fromtimestamp(imTime).strftime('%Y-%m-%d-%H-%M-%S-%f'))
#            _ = gettext.gettext                                                                                              
            for i,cnt3 in enumerate([imageDir]):
#                logger.info(_('File %i'), i)
                if overWrite or not os.path.isfile(os.path.normpath(os.path.join(cnt3,imFilename))):
                    fid.seek(offset,0)
                    im = np.fromfile(fid,dtype='<u1',count=imageLengthBytes).reshape(imageSizeNy,imageSizeNx)
                    if len(cropcoords)!=0:
                        im = im[cropcoords[0]:cropcoords[1]][cropcoords[2]: cropcoords[3]]
                    print("{}, {}\n".format(cnt2,imFilename))
                    PIL.Image.fromarray(im).save(os.path.normpath(os.path.join(cnt3,imFilename)),optimize=True)
                else:
                    continue

#        logger.info(_('Done!'))
    # add the seq file to the blacklist in the image dir to mark that it has already been done
    _ , seqFilename = os.path.split(thisSequenceFilePath)
    addToBlacklist(imageDir, seqFilename)
    
def loadBlacklist(thisdir):
    if not os.path.isdir(thisdir):
        print("Cannot load a blacklist as {} is not a directory".format(thisdir))
    if not os.path.isfile(os.path.normpath(os.path.join(thisdir,'blacklist.txt'))):
        print("No {} file found.".format(os.path.normpath(os.path.join(thisdir,'blacklist.txt'))))
        with open(os.path.normpath(os.path.join(thisdir,'blacklist.txt')),'w+') as fid:
            pass
    with open(os.path.normpath(os.path.join(thisdir,'blacklist.txt')),'r') as fid:
        text = fid.read().splitlines()
    
    return text

def addToBlacklist(thisdir,files):
    if not os.path.isdir(thisdir):
        print("Cannot save to blacklist as {} is not a directory".format(thisdir))
    if isinstance(files,str):
        files = [files]
    existFiles = loadBlacklist(thisdir)
    files = [x for x in files if x not in existFiles]
    if not files: return
    while not open(os.path. normpath(os.path.join(thisdir,'blacklist.txt')),'a'):
        pass                                
    with open(os.path.normpath(os.path.join(thisdir,'blacklist.txt')),'a') as fid:
        for file in files:
            fid.write("\n".join(file))

def main():
    root = tk.Tk()
    root.withdraw()
#    logging.basicConfig()
#    logging.getLogger().setLevel(logging.INFO)
    NorPixSequenceFileImageConverter()
    
if __name__ == '__main__':
    main()
