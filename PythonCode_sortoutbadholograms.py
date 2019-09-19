'''
Read in holograms and find bad ones
'''
import os
import shutil
import numpy as np
import datetime
import re
import matplotlib as mpl
from matplotlib.backends.backend_tkagg  import  FigureCanvasTkAgg
import tkinter as tk
#import sys
import PIL

class Sortoutbadholograms(object):
    '''This class is ment to help sorting out bad holograms'''
    def __init__(self, holodir):
        self.badind = []
        self.holofiles = []
        self.holobyts = []
        self.holomi = []
        self.holosi = []
        self.holomx = []
        self.holosx = []
        self.holomy = []
        self.holosy = []
        self.badholofiles = []
        self.unreadableholofiles = []
        self.holotimes = []
        self._holonos = None

        # Get the filenames
        tholofiles = [x for x in os.scandir(holodir) if x.name.endswith('.png')]
        tholobyts = [x.stat().st_size for x in tholofiles]
        tholofiles = np.array([x.path for x in tholofiles])

        # Get statistics from holograms and filtered holograms
        tbadholofiles = np.array([False]*np.size(tholofiles))
        tunreadableholofiles = np.array([False]*np.size(tholofiles))

        # Collect any precalculated statistics left from
        # NorpixSequenceFileConverter.
        tholofiles, existingstats = self.getexistingstats(holodir, tholofiles)
        if not existingstats:
#            existingstats = np.empty(shape=(len(tholofiles), 6), dtype='f8')
#            existingstats = [None]*len(tholofiles)
#            for nr, name in enumerate(tholofiles):
#                existingstats[nr] = self.basicstats(name)
            existingstats = [self.basicstats(x) for x in tholofiles]
        existingstats = np.asanyarray(existingstats)
        # Determine what's left to calculate
        dooverinds = np.argwhere(np.any(np.isnan(np.array(existingstats.tolist())), axis=1))
        #uS = etd(clock, 1, numel(dooverinds), 60)
        for ind in dooverinds:
            #uS = etd(uS, cnt)
            #starttime=clock
#            try:
#                with open(tholofiles[ind], 'rb') as fid:
#                    im =np.fromfile(fid, dtype='f4')
#            except:
#                tunreadableholofiles[ind] = True
#                continue
            existingstats[ind] = self.basicstats(tholofiles[ind])
             #fprintf('%d in %.1f seconds\n',ind, etime(clock, starttime))

        # save the calculated variables
        self.holofiles = tholofiles
        self.holobyts = tholobyts
        [setattr(self, key, existingstats[key]) for key in existingstats.dtype.names]
        self.badholofiles = tbadholofiles
        self.unreadableholofiles = tholofiles[tunreadableholofiles]

        self.markbaddies()

#        save([holodir filesep 'sortoutbaddies_holostats.mat'], 'self')

    @property
    def holonos(self):
        return range(len(self.holofiles))

    def plotmisi(self):
        root = tk.Tk()
        MyWindow(root, self).pack()
        root.mainloop()
#        fig156=plt.figure(num=156, clear=True)
#        contourPlot(self.holosi/self.holomi, self.holomi,'PrcTile',[0 100])
#        plt.hold('on')
#        scatter(self.holosi(self.badind)/self.holomi(self.badind), self.holomi(self.badind),'rp')
#        hold off
#        ylabel('Intensity Mean')
#        xlabel('Intensity Dispersion/Intensity Mean')


    def plotmxmy(self):
        root = tk.Tk()
        MyWindow(root, self).pack()
        root.mainloop()
#        fig157=plt.figure(num=157, clear=True)
#        contourPlot(self.holomx, self.holomy, 'PrcTile',[0 100])
#        hold on
#        scatter(self.holomx(self.badind), self.holomy(self.badind),'rp')
#        hold off
#        xlabel('X-Position Mean')
#        ylabel('Y-Position Mean')


    def plotsxsy(self):
        root = tk.Tk()
        MyWindow(root).pack()
        root.mainloop()
#        fig158=plt.figure(num=158, clear=True)
#        contourPlot(self.holosx, self.holosy, 'PrcTile',[0 100])
#        hold on
#        scatter(self.holosx(self.badind), self.holosy(self.badind),'rp')
#        xlabel('X-Position Dispersion')
#        ylabel('Y-Position Dispersion')
#        hold off


    def plotmisiseries(self):
        root = tk.Tk()
        MyWindow(root).pack()
        root.mainloop()
#        tholotimes = self.getholotimes
#        fig2000=plt.figure(num=2009, clear=True)
#        subplot(2, 1, 1)
#        plot(tholotimes, self.holosi,'g.',tholotimes(self.badind),self.holosi(self.badind),'b.')
#        title('Dispersion in intensity')
#        ylabel('Intensity')
#        if self.isholotime(tholotimes): plt.datetick('x')
#        subplot(2, 1, 2)
#        plot(tholotimes, self.holomi,'g.',tholotimes(self.badind),self.holomi(self.badind),'b.')
#        title('Mean in intensity')
#        if self.isholotime(tholotimes):
#            datetick('x') xlabel('Hologram Time')
#        else:
#            xlabel('Hologram number')
#        ylabel('Intensity')
#        figure(gcf)


    def plotmxsxseries(self):
        root = tk.Tk()
        MyWindow(root).pack()
        root.mainloop()
#        tholotimes = self.getholotimes
#        fig2010=plt.figure(2010)
#        subplot(2, 1, 1)
#        plot(tholotimes, self.holosx,'g.',tholotimes(self.badind),self.holosx(self.badind),'b.')
#        title('Dispersion in x direction')
#        ylabel('Intensity')
#        if self.isholotime(tholotimes), datetick('x')
#        subplot(2, 1, 2)
#        plot(tholotimes, self.holomx,'g.',tholotimes(self.badind),self.holomx(self.badind),'b.')
#        title('Mean in x position')
#        ylabel('Intensity')
#        if self.isholotime(tholotimes), datetick('x') xlabel('Hologram Time')
#        else, xlabel('Hologram number')
#        figure(gcf)


    def plotmysyseries(self):
        root = tk.Tk()
        MyWindow(root).pack()
        root.mainloop()
#        tholotimes = self.getholotimes
#        fig2011=plt.figure(2011)
#        subplot(2, 1, 1)
#        plot(tholotimes, self.holosy,'g.',tholotimes(self.badind),self.holosy(self.badind),'b.')
#        title('Dispersion in y direction')
#        ylabel('Intensity')
#        if self.isholotime(tholotimes), datetick('x')
#        subplot(2, 1, 2)
#        plot(tholotimes, self.holomy,'g.',tholotimes(self.badind),self.holomy(self.badind),'b.')
#        title('Mean in y position')
#        ylabel('Intensity')
#        if self.isholotime(tholotimes), datetick('x') xlabel('Hologram Time')
#        else, xlabel('Hologram number')
#        figure(gcf)


    def getholotimes(self):
        if self.holotimes:
            timestamps = self.holotimes
            return
        def maketimestamp(filename):
            # Make the timestamp from the filename
            temp = re.findall('_([0-9-]+)', filename)
            tempdate = datetime.datetime.fromordinal(temp, 'yyyy-mm-dd-HH-MM-SS-FFF')
            return datetime.datetime.fromtimestamp(tempdate)
        try:
            timestamps = [maketimestamp(x) for x in self.holofiles]
            self.holotimes = timestamps
        except:
            print('Filenames of holograms not suitable for holotimes conversion!')
            timestamps = self.holonos


    def markbaddies(self):
        '''determine mean intensity and set limit to 5% above and below'''
        mi = np.mean(self.holomi)
        minMI = mi*(1-0.05)
        maxMI = mi*(1+0.05)

        self.badind = (self.holomi < minMI) | (self.holomi > maxMI)
        self.badholofiles = self.holofiles[self.badind]


    def movebaddies(self, inselfdir):
        '''move bad holograms found by sortoutbadholograms'''
        baddies = np.vstack((self.badholofiles, self.unreadableholofiles))
        # if folder does not exist create it
        if not os.path.isdir(inselfdir):
            os.mkdir(inselfdir)
        # move the bad holograms
        #uS = etd(clock, 1, length(baddies), 60)
        for i_i in range(len(baddies)):
            if os.path.isfile(baddies[i_i]):
                shutil.move(baddies[i_i], inselfdir)
            #uS = etd(uS, i_i)

    @staticmethod
    def isholotime(nums):
        '''Determin if nums is a valid holotime'''
        return not np.all(nums == np.fix(nums))

    @staticmethod
    def basicstats(im):
        '''Calculate descriptive statistics for image im'''
        if isinstance(im, str):
            with open(im, 'rb') as fid:
                im = np.asarray(PIL.Image.open(fid))
        mi = np.mean(im)
        si = np.std(im)
        Ny, Nx = im.shape
        sum_ = np.sum(im)
        #mx = sum(sum(bsxfun(@times, 1:Nx, im)))/sum(im(:))
        #my = sum(sum(bsxfun(@times,(1:Ny)',im)))/sum(im(:))
        mx = np.sum(np.sum(np.array((range(Nx))*im)))/sum_
        my = np.sum(np.sum(np.array((range(Ny))*im.T)))/sum_
        #sx = np.sqrt(sum(sum(bsxfun(@times,((1:Nx)-mx).^2, im)))/sum(im(:)))
        #sy = np.sqrt(sum(sum(bsxfun(@times,((1:Ny)-my).^2, im)))/sum(im(:)))
        sx = np.sqrt(np.sum(np.sum((np.array((range(Nx)))-mx)**2 * im))/sum_)
        sy = np.sqrt(np.sum(np.sum((np.array((range(Ny)))-my)**2 * im.T))/sum_)
        names = ['holomi', 'holosi', 'holomx', 'holomy', 'holosx', 'holosy']
        formats = ['f8', 'f8', 'f8', 'f8', 'f8', 'f8']
        ddtype = dict(names=names, formats=formats)
        return np.array((mi, si, mx, my, sx, sy), dtype=ddtype)
#        return np.array([mi, si, mx, my, sx, sy])
#        return {'mi':mi, 'si':si, 'mx':mx, 'my':my, 'sx':sx, 'sy':sy}

    @staticmethod
    def getexistingstats(holodir, tholofiles):
        '''Search for and load precalculated statistics from the holograms done
         from NorpixSequenceFileConverter.'''
        # Get the filenames of any precalculated stats files.
        statfiles = [x.path for x in os.scandir(holodir) if x.name.startswith(
            'sortoutthebaddies_') and x.name.endswith('_holostats.mat')]
        # Load the precalculated stat files
        jojo = np.array([None]*len(statfiles))
        if jojo.size == 0:
            return tholofiles, None
        for cnt, file in enumerate(statfiles):
            with open(file, 'rb') as fid:
                jojo[cnt] = fid.load()
        # Collapse the data to a single struct
        #jojo = [jojo{:}] jojo = [jojo.meme]  #fpp
        jojoholofiles = [x for x in jojo.holofilenames]
        #jojoholofiles = vertcat(jojoholofiles{:})
        jojoholomi = np.array(jojo.holomi)
        jojoholosi = np.array(jojo.holosi)
        jojoholomx = np.array(jojo.holomx)
        jojoholomy = np.array(jojo.holomy)
        jojoholosx = np.array(jojo.holosx)
        jojoholosy = np.array(jojo.holosy)
        # Make sure the values are sorted
        #jojoholofiles = np.unique(jojoholofiles)
        #jojoholomi = np.take(jojoholomi, inds)
        #jojoholosi = np.take(jojoholosi, inds)
        #jojoholomx = np.take(jojoholomx, inds)
        #jojoholomy = np.take(jojoholomy, inds)
        #jojoholosx = np.take(jojoholosx, inds)
        #jojoholosy = np.take(jojoholosy, inds)

        _, tholofiles = [os.path.split(x) for x in tholofiles]
        _, jojoholofiles = [os.path.split(x) for x in jojoholofiles]
        # Find where to insert these values in the whole tholofiles ro
        #tholofiles, inds = np.union1d((tholofiles, jojoholofiles),return_index=True)
        #tholomi = np.nan(np.size(tholofiles))
        #tholosi = np.nan(np.size(tholofiles))
        #tholomx = np.nan(np.size(tholofiles))
        #tholomy = np.nan(np.size(tholofiles))
        #tholosx = np.nan(np.size(tholofiles))
        #tholosy = np.nan(np.size(tholofiles))
        tholomi = jojoholomi
        tholosi = jojoholosi
        tholomx = jojoholomx
        tholomy = jojoholomy
        tholosx = jojoholosx
        tholosy = jojoholosy
        return tholofiles, {'tholomi' : tholomi, 'tholes' : tholosi, 'tholomx' : tholomx,
                            'tholomy' : tholomy, 'tholosx' : tholosx, 'tholosy' : tholosy}

class MyWindow(tk.Frame):
    def __init__(self, parent, cls):
        tk.Frame.__init__(self, parent)
        self.fig = mpl.figure.Figure((6, 6), dpi=100)
        canvas = FigureCanvasTkAgg(self.fig, master=self)
        canvas.get_tk_widget().grid(row=0, column=0, rowspan=2)
        label = tk.Label(self, text="Choose Type")
        label.grid(row=0, column=1)
        choosebox = tk.Listbox(self)
        [choosebox.insert(x, str(y)) for x, y in enumerate(
                ['plotmisi', 'plotmxmy','plotsxsy', 'plotmisiseries',
                 'plotmxsxseries', 'plotmysyseries'])]
#        choosebox.bind('<<ListboxSelect>>', getattr(cls,choosebox.get(ACTIVE)())
        choosebox.grid(row=1, column=1)
        print(choosebox.get(active))

def main():
    """Test doc String from NorpixSequenceFileConverter"""
#    root = tk.Tk()
#    root.withdraw()
    test = Sortoutbadholograms(r'C:\Users\mbexwws2\Documents\HOLOSUITE\holograms')
    test.getholotimes()
if __name__ == '__main__':
    main()

