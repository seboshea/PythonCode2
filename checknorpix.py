# -*- coding: utf-8 -*-
"""
Small program to control if norpixconverter did well. Should be impolemented 
into norpixconverter, maybe as additional information in the blacklist file
to avoid multiple runs of the fils which takes quite long.

Created on Thu Apr  5 07:29:02 2018

@author: mbexwws2
"""
import os
import tkinter.filedialog as fd
import multiprocessing
import functools
import datetime
import tkinter as tk
import tqdm
import numpy as np


def checknorpix(seqpath='', imagedir='', imageprefix='C075_'):
    
    root = tk.Tk()
    root.withdraw()
    
    if not seqpath or not os.path.isdir(seqpath):
        seqdir = fd.askdirectory(title='Please choose a folder for seq files')
        if not seqdir:
            tqdm.print('User failed to choose a valid folder for *.seq files\n \
                  Operation cancelled by user')
            return
    else:
        seqdir = seqpath
    seqpaths = [file.path for file in os.scandir(os.path.normpath(os.path.join(seqdir)))
                if file.path.endswith('.seq')]
    if not imagedir or not os.path.isdir(imagedir):
        imagedir = fd.askdirectory(title='Please choose a flolder for holograms')
        if not imagedir:
            print('User failed to choose a valid folder for holograms \
                  production\n Operation cancelled by user')
            return
    pngpaths = [file.name for file in os.scandir(os.path.normpath(os.path.join(imagedir)))
                if file.path.endswith('.png')]

    isfinished = False

    while not isfinished:
        isfinished = proc_seqfiles(seqdir, seqpaths, imagedir, pngpaths, imageprefix)
        
def proc_seqfiles(seqdir, seqpaths, imagedir, pngpaths, imageprefix):
    """Function to process all Seq files in folder"""
    for cnt, file in enumerate(tqdm.tqdm(seqpaths)):
        process_oneseqfile(file, imagedir, pngpaths, imageprefix)
#    with multiprocessing.Pool() as pool:
#        func = functools.partial(process_oneseqfile, pngpaths=pngpaths,
#                                 imageprefix=imageprefix)
#        pool.map(func, seqpaths)

    return True

def process_oneseqfile(thisseqpath, imagedir, pngpaths, imageprefix):
    """
    Helper function of proc_seqfiles to process individual Seq files
    """
    startt = datetime.datetime.now()
    filesize = os.path.getsize(thisseqpath)
    successful = 0
    
    with open(thisseqpath, 'rb') as fid:
        temp = np.fromfile(fid, dtype='<u4', count=155, sep='')
        magicnumber = temp[0]
        if magicnumber != 65261:
            tqdm.tqdm.write("{} has the wrong magic number for a sequence file".format(thisseqpath))
            return None
        fileheaderlength = 8192

        imageSizeNx, imageSizeNy = temp[137:139]
        imageBitDepth, imageBitDepthReal = temp[139:141]
        imageLengthBytes, imageFormat = temp[141:143]

        # Get the number of images that are "supposed" to be in this file
        numImagesShould = temp[143]
        if (imageBitDepth != 8) | (imageBitDepthReal != 8) | (imageFormat != 100):
            tqdm.tqdm.write("This function only supports Monochrome 8-bit images \n \
                            Error in {}".format(thisseqpath))
            endt = datetime.datetime.now()-startt
            towrite = "\t".join([os.path.basename(thisseqpath), str(1),
                         str(np.NaN), str(np.NaN), str(np.NaN), str(startt), str(endt),
                         str(0), ",".join([])])
    
            addToChecklist(imagedir, "{}{}".format(imageprefix, '_holograms_check.txt'),
                           towrite)
            return None
        # Get the spacing between images as they are separated by the image size
        # + the image footer up to the next sector boundary
        imageSpacing = temp[145]
        
        shouldBeNumImages = int((filesize - fileheaderlength)/(imageSpacing))
        if numImagesShould != shouldBeNumImages:
            numImages = shouldBeNumImages
        else:
            numImages = numImagesShould
        
        offset = fileheaderlength + range(numImages)*imageSpacing
        try:
            offset_diff = np.append(offset[0]+imageLengthBytes, np.diff(offset))-imageLengthBytes
        except:
            successful = 1
        Time1 = np.empty((numImages,), dtype='<u4')
        Time2 = np.empty((numImages,), dtype='<u2')
        Time3 = np.empty((numImages,), dtype='<u2')
        #imFilename = np.array(['']*(numImages), dtype=object)
        ddtype = np.dtype([('f0', '<u4'), ('f1', '<u2'), ('f2', '<u2')])

        fid.seek(8, 0)
        for cnt2 in range(numImages):
            o = offset_diff[cnt2]
            fid.seek(o+imageLengthBytes-8, 1)
            Time1[cnt2], Time2[cnt2], Time3[cnt2] = np.fromfile(fid, dtype=ddtype,
                                                                count=1, sep='')[0]
        imTime = Time1 + Time2/1000 + Time3/1e6
        imFilename = np.array(["{}{}.png".format(
                imageprefix, datetime.datetime.fromtimestamp(x).strftime('%Y-%m-%d-%H-%M-%S-%f'))
                               for x in imTime if x], dtype=object)
        notexistingfiles = [file for file in imFilename if file not in pngpaths]
    # add the seq file to the blacklist in the image dir to mark that it has 
    # already been done
    endt = datetime.datetime.now()-startt
    
    towrite = "\t".join([os.path.basename(thisseqpath), str(successful),
                         str(numImagesShould), str(shouldBeNumImages),
                         str(len(imFilename)), str(startt), str(endt),
                         str(len(notexistingfiles)), ",".join(notexistingfiles)])
    
    addToChecklist(imagedir, "{}{}".format(imageprefix, '_holograms_check.txt'),
                   towrite)

def loadChecklist(thisdir):
    """Function to load a blacklist of already used Seq file names"""
    if not os.path.isdir(os.path.dirname(thisdir)):
        print("Cannot load a checklist as {} is not a directory".format(thisdir))
    if not os.path.isfile(thisdir):
        print("No {} file found.".format(thisdir))
        with open(thisdir, 'w+') as fid:
            header = "\t".join(['SeqFile', 'Error', '# reported', '# should (size)',
                               '# is', 'starttime', 'endtime',
                               '# missing', 'Png Names'])
            fid.write(header+'\n')
    with open(thisdir, 'r') as fid:
        header = fid.readline()
        text = fid.read().splitlines()
    lines = [l.split('\t') for l in text]
    if not lines:
        return ['']
    text, *reason = zip(*lines)

    return text, reason

def addToChecklist(thisdir='.', checkname='Checklist.txt', inputtxt=''):
    """Function to add name of Seq file to a blacklist"""
    checkpath = os.path.normpath(os.path.join(thisdir, checkname))
    if not os.path.isdir(thisdir):
        print("Cannot save to checklist as {} is not in directory".format(checkpath))
    existFiles, *_ = loadChecklist(checkpath)
    while not open(checkpath, 'a'):
        pass
    with open(checkpath, 'a') as fid:
        fid.write(inputtxt+'\n')
        
def main():
    """Test doc String from NorpixSequenceFileConverter"""
    root = tk.Tk()
    root.withdraw()
#    checknorpix(imageprefix='C164_')
    checknorpix(seqpath='K:\\C161\\HALOHOLO\\',imagedir='K:\\C161\\HALOHOLO\\',imageprefix='C161_')

if __name__ == '__main__':
    __spec__=None
    main()
