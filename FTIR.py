#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from PyQt5.QtWidgets import (QWidget, QPushButton,QLabel,QFileDialog,QComboBox,
    QHBoxLayout, QVBoxLayout, QGridLayout, QLineEdit, QApplication,QRadioButton,QErrorMessage)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from scipy.optimize import minimize


class FTIRAnalysis(QWidget):

    sampledata=None
    refdata=None
    diffFCA=None
    diffOLine=None
    diffCLine=None
    poptFCA=None
    sample_thick=None
    PurespectrumFCA=None
    PureOxygen=None
    PureCarbon=None
    refdatasubC=None
    sampledatasubC=None
    refdatasubO=None
    sampledatasubO=None
    thicknessfactorLine=None
    thicknessfactorFCA=None
    thicknessfactor=None


    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):


        SampleSection = QVBoxLayout()

        SampleSection1= QHBoxLayout()
        LoadSampleButton = QPushButton('Load Sample')
        LoadSampleButton.clicked.connect(self.loadsample)
        LSampleThick = QLabel('Thickness [\u03BCm]')
        LSample = QLabel('Sample')
        SampleSection1.addStretch(1)
        SampleSection1.addWidget(LSample)
        SampleSection1.addWidget(LoadSampleButton)
        SampleSection1.addStretch(1)
        SampleSection1.addWidget(LSampleThick)

        SampleSection2= QHBoxLayout()
        self.SampleName = QLineEdit()
        self.SampleName.setMinimumWidth(300)
        self.SampleThick = QLineEdit()
        self.SampleThick.textChanged[str].connect(self.ThickChange)
        SampleSection2.addWidget(self.SampleName)
        SampleSection2.addWidget(self.SampleThick)

        SampleSection.addLayout(SampleSection1)
        SampleSection.addLayout(SampleSection2)

        RefSection = QVBoxLayout()

        RefSection1= QHBoxLayout()
        LRef = QLabel('Reference')
        LoadRefButton = QPushButton('Load Reference')
        LoadRefButton.clicked.connect(self.loadref)
        LRefThick = QLabel('Thickness [\u03BCm]')
        RefSection1.addStretch(1)
        RefSection1.addWidget(LRef)
        RefSection1.addWidget(LoadRefButton)
        RefSection1.addStretch(1)
        RefSection1.addWidget(LRefThick)

        RefSection2= QHBoxLayout()
        self.RefName = QLineEdit()
        self.RefName.setMinimumWidth(300)
        self.RefThick = QLineEdit()
        self.RefThick.textChanged[str].connect(self.ThickChange)
        RefSection2.addWidget(self.RefName)
        RefSection2.addWidget(self.RefThick)

        RefSection.addLayout(RefSection1)
        RefSection.addLayout(RefSection2)

        DataSection=QHBoxLayout()
        DataSection.addLayout(SampleSection)
        DataSection.addLayout(RefSection)

        PhononCorrectSection = QHBoxLayout()

        ThickSection = QVBoxLayout()
        Lthickfactor = QLabel('Thickness Ratio')
        self.Thickfactor = QLineEdit()
        ThickSection.addWidget(Lthickfactor)
        ThickSection.addWidget(self.Thickfactor)
        self.Thickfactor.setReadOnly(True)

        CorrectSection = QGridLayout()
        self.ChooseSL = QRadioButton('Straght lines Correction')
        self.ChooseSL.setChecked(True)
        self.ChooseSL.toggled.connect(self.Choosechange)
        self.ChooseFCA = QRadioButton('FCA Correction')
        self.GetSLfactor = QPushButton('Get Phonon factor from Straight lines correction')
        self.GetSLfactor.setDisabled(True)
        self.GetSLfactor.clicked.connect(self.getfactorlineButtom)
        self.GetFCAfactor = QPushButton('Get Phonon factor from FCA correction')
        self.GetFCAfactor.setDisabled(True)
        self.GetFCAfactor.clicked.connect(self.getfactorFCAButtom)
        self.SLfactor = QLineEdit()
        self.SLfactor.textChanged[str].connect(self.FactorChange)
        self.FCAfactor = QLineEdit()
        self.FCAfactor.textChanged[str].connect(self.FactorChange)
        CorrectSection.addWidget(self.ChooseSL,0,0)
        CorrectSection.addWidget(self.ChooseFCA,0,1)
        CorrectSection.addWidget(self.GetSLfactor,1,0)
        CorrectSection.addWidget(self.GetFCAfactor,1,1)
        CorrectSection.addWidget(self.SLfactor,2,0)
        CorrectSection.addWidget(self.FCAfactor,2,1)

        PhononCorrectSection.addLayout(ThickSection)
        PhononCorrectSection.addLayout(CorrectSection)

        CalcuSection = QVBoxLayout()
        self.CalOxygen = QPushButton('Calculate Oxygen')
        self.CalOxygen.clicked.connect(self.CalculateOxygen)
        self.CalOxygen.setDisabled(True)
        self.CalCarbon = QPushButton('Calculate Carbon')
        self.CalCarbon.setDisabled(True)
        self.CalCarbon.clicked.connect(self.CalculateCarbon)
        CalcuSection.addWidget(self.CalOxygen)
        CalcuSection.addWidget(self.CalCarbon)

        h1box = QHBoxLayout()
        h1box.addLayout(DataSection)
        h1box.addStretch(1)
        h1box.addLayout(PhononCorrectSection)
        h1box.addLayout(CalcuSection)

        SpectrumPlot=QVBoxLayout()
        self.figure1 = plt.figure()
        self.canvas1 = FigureCanvas(self.figure1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        SpectrumPlot.addWidget(self.canvas1)
        SpectrumPlot.addWidget(self.toolbar1)
        self.ax1 = self.figure1.add_subplot(111)
        self.ax1.set_xlabel(r'Wavenumber $[cm^{-1}]$')
        self.ax1.set_ylabel('Absorbance')
        self.ax1.set_xlim([1800,400])
        self.figure1.tight_layout()

        ZoominPlot=QVBoxLayout()
        self.figure2 = plt.figure()
        self.canvas2 = FigureCanvas(self.figure2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        ZoominPlot.addWidget(self.canvas2)
        ZoominPlot.addWidget(self.toolbar2)
        self.ax2 = self.figure2.add_subplot(211)
        self.ax2.set_title('Oxygen')
        self.ax2.set_xlim([1250,1000])
        self.ax3 = self.figure2.add_subplot(212)
        self.ax3.set_xlabel(r'Wavenumber $[cm^{-1}]$')
        self.ax3.set_title('Carbon')
        self.ax3.set_xlim([700,500])
        self.figure2.tight_layout()

        h2box = QHBoxLayout()
        h2box.addLayout(SpectrumPlot)
        h2box.addLayout(ZoominPlot)

        vbox = QVBoxLayout()
        vbox.addLayout(h1box)
        vbox.addLayout(h2box)

        self.setLayout(vbox)

        self.setGeometry(300, 200, 1400, 700)
        self.setWindowTitle('FTIR for Oxygen and Carbon concentrations')
        self.show()

    def loadsample(self):
        Samplename = QFileDialog.getOpenFileName(self, caption='Choose Sample file',filter='*.csv')
        if Samplename[0] != '':
            self.sampledata=np.genfromtxt(Samplename[0],delimiter=',')
            if self.sampledata.shape[1] != 2:
                self.sampledata = None
                self.SampleName.setText('Invalid Data File!')
                self.SampleName.setStyleSheet("color:red")
            elif (self.refdata is not None) and (self.sampledata.shape[0] != self.refdata.shape[0] +1):
                self.sampledata = None
                self.diffFCA = None
                self.poptFCA = None
                self.PurespectrumFCA = None
                self.diffCLine = None
                self.diffOLine = None
                self.refdata = None
                self.refdatasubC = None
                self.refdatasubO = None
                self.sampledatasubC = None
                self.sampledatasubO = None
                self.RefName.clear()
                self.RefThick.clear()
                self.SampleName.clear()
                self.SampleThick.clear()
                self.plotdata()
                error_dialog = QErrorMessage()
                error_dialog.showMessage('The length of sample and reference spectrum does not match!')
                error_dialog.exec_()
            else:
                self.SampleThick.setText((os.path.basename(Samplename[0])).split()[0])
                self.sampledata=self.sampledata[1:,:]
                index = np.argsort(self.sampledata[:,0])
                self.sampledata=self.sampledata[index,:]
                self.SampleName.setStyleSheet("color:black")
                self.SampleName.setText(os.path.basename(Samplename[0]))
                self.diffFCA = None
                self.poptFCA = None
                self.PurespectrumFCA = None
                self.diffCLine = None
                self.diffOLine = None
                self.FCAfactor.clear()
                self.SLfactor.clear()
                self.subtractlineSample()
                self.plotdata()
                if self.refdata is not None:
                    self.GetSLfactor.setDisabled(False)
                    self.GetFCAfactor.setDisabled(False)
                    if self.sample_thick is not None:
                        self.CalOxygen.setDisabled(False)
                        self.CalCarbon.setDisabled(False)

    def loadref(self):
        Refname = QFileDialog.getOpenFileName(self, caption='Choose Reference file',filter='*.csv')
        if Refname[0] != '':
            self.refdata=np.genfromtxt(Refname[0],delimiter=',')
            if self.refdata.shape[1] != 2:
                self.refdata = None
                self.RefName.setText('Invalid Data File!')
                self.RefName.setStyleSheet("color:red")
            elif (self.sampledata is not None) and (self.sampledata.shape[0] != self.refdata.shape[0]-1):
                self.sampledata = None
                self.diffFCA = None
                self.poptFCA = None
                self.PurespectrumFCA = None
                self.diffCLine = None
                self.diffOLine = None
                self.refdata = None
                self.refdatasubC = None
                self.refdatasubO = None
                self.sampledatasubC = None
                self.sampledatasubO = None
                self.RefName.clear()
                self.RefThick.clear()
                self.SampleName.clear()
                self.SampleThick.clear()
                self.plotdata()
                error_dialog = QErrorMessage()
                error_dialog.showMessage('The length of sample and reference spectrum does not match!')
                error_dialog.exec_()
            else:
                self.RefThick.setText((os.path.basename(Refname[0])).split()[0])
                self.refdata=self.refdata[1:,:]
                index = np.argsort(self.refdata[:,0])
                self.refdata=self.refdata[index,:]
                self.RefName.setStyleSheet("color:black")
                self.RefName.setText(os.path.basename(Refname[0]))
                self.diffFCA = None
                self.poptFCA = None
                self.PurespectrumFCA = None
                self.diffCLine = None
                self.diffOLine = None
                self.FCAfactor.clear()
                self.SLfactor.clear()
                self.subtractlineRef()
                self.plotdata()
                if self.refdata is not None:
                    self.GetSLfactor.setDisabled(False)
                    self.GetFCAfactor.setDisabled(False)
                    if self.sample_thick is not None:
                        self.CalOxygen.setDisabled(False)
                        self.CalCarbon.setDisabled(False)

    def getfactorlineButtom(self):
        self.FitThicknessfactorL()

    def getfactorFCAButtom(self):
        self.FitThicknessfactorFCA()

    def Choosechange(self):
        self.plotdata()

    def ThickChange(self):
        try:
            self.sample_thick = float(self.SampleThick.text())
            if (self.refdata is not None) and (self.sampledata is not None):
                self.CalOxygen.setDisabled(False)
                self.CalCarbon.setDisabled(False)
        except ValueError:
            self.sample_thick = None
            self.thicknessfactor = None
            self.CalOxygen.setDisabled(True)
            self.CalCarbon.setDisabled(True)
            self.Thickfactor.clear()
            self.SampleThick.clear()
        try:
            ref_thick = float(self.RefThick.text())
            if self.sample_thick is not None:
                self.Thickfactor.setText(str(round(self.sample_thick/ref_thick,3)))
                self.thicknessfactor = self.sample_thick/ref_thick
        except ValueError:
            self.thicknessfactor = None
            self.Thickfactor.clear()
            self.RefThick.clear()


    def FactorChange(self):
        try:
            self.thicknessfactorLine = float(self.SLfactor.text())
            self.substractreferenceline()
        except ValueError:
            self.thicknessfactorLine = None
            self.SLfactor.clear()
        try:
            self.thicknessfactorFCA = float(self.FCAfactor.text())
            self.substractreferenceFCA()
            self.fitFCA()
            self.substractFCA()
        except ValueError:
            self.thicknessfactorFCA = None
            self.FCAfactor.clear()
        self.plotdata()

    def subtractlineRef(self):
        index550 = (np.abs(self.refdata[:,0]-550)).argmin()
        index650 = (np.abs(self.refdata[:,0]-650)).argmin()
        index1020 = (np.abs(self.refdata[:,0]-1020)).argmin()
        index1220 = (np.abs(self.refdata[:,0]-1220)).argmin()
        Ref550 = np.average(self.refdata[index550-3:index550+3,1])
        Ref650 = np.average(self.refdata[index650-3:index650+3,1])
        refliene = self.refdata[:,0]*(Ref650-Ref550)/100+55/10*(Ref550*65/55-Ref650)
        self.refdatasubC = self.refdata[:,1]-refliene
        Ref1020 = np.average(self.refdata[index1020-3:index1020+3,1])
        Ref1220 = np.average(self.refdata[index1220-3:index1220+3,1])
        refliene = self.refdata[:,0]*(Ref1220-Ref1020)/200+1020/200*(Ref1020*1220/1020-Ref1220)
        self.refdatasubO = self.refdata[:,1]-refliene

    def subtractlineSample(self):
        index550 = (np.abs(self.sampledata[:,0]-550)).argmin()
        index650 = (np.abs(self.sampledata[:,0]-650)).argmin()
        index1020 = (np.abs(self.sampledata[:,0]-1020)).argmin()
        index1220 = (np.abs(self.sampledata[:,0]-1220)).argmin()
        sample550 = np.average(self.sampledata[index550-3:index550+3,1])
        sample650 = np.average(self.sampledata[index650-3:index650+3,1])
        sampleliene = self.sampledata[:,0]*(sample650-sample550)/100+55/10*(sample550*65/55-sample650)
        self.sampledatasubC = self.sampledata[:,1]-sampleliene
        sample1020 = np.average(self.sampledata[index1020-3:index1020+3,1])
        sample1220 = np.average(self.sampledata[index1220-3:index1220+3,1])
        sampleliene = self.sampledata[:,0]*(sample1220-sample1020)/200+1020/200*(sample1020*1220/1020-sample1220)
        self.sampledatasubO = self.sampledata[:,1]-sampleliene

    def substractreferenceline(self):
        if (self.sampledata is not None) and (self.refdata is not None) and (self.thicknessfactorLine is not None):
            self.diffCLine=self.sampledatasubC-self.refdatasubC*self.thicknessfactorLine
            self.diffOLine=self.sampledatasubO-self.refdatasubO*self.thicknessfactorLine

    def FitThicknessfactorL(self):
        try:
            index615=(np.abs(self.sampledata[:,0]-615)).argmin()
            index621=(np.abs(self.sampledata[:,0]-621)).argmin()
            self.thicknessfactorLine=np.sum(self.sampledatasubC[index615:index621])/np.sum(self.refdatasubC[index615:index621])
            self.SLfactor.setText(str(self.thicknessfactorLine))
        except:
            pass

    def substractreferenceFCA(self):
        if (self.sampledata is not None) and (self.refdata is not None) and (self.thicknessfactorFCA is not None):
            self.diffFCA=self.sampledata[:,1]-self.refdata[:,1]*self.thicknessfactorFCA

    def FitThicknessfactorFCA(self):
        try:
            index1 = self.refdata[:,0]>500
            index1 *= self.refdata[:,0]<590
            index2 = self.refdata[:,0]>610
            index2 *= self.refdata[:,0]<1020
            index3 = self.refdata[:,0]>1200
            index3 *= self.refdata[:,0]<1600
            index = index1+index2+index3
            def residual(k):
                diff = self.sampledata[index,1]-self.refdata[index,1]*k
                p = np.polyfit(1/self.refdata[index,0],diff,2)
                return np.sum(abs(diff-np.polyval(p,1/self.refdata[index,0])))
            res = minimize(residual,1,method='Nelder-Mead', tol=1e-6)
            self.thicknessfactorFCA = res.x[0]
            self.FCAfactor.setText(str(self.thicknessfactorFCA))
        except:
            pass

    def fitFCA(self):
        if (self.refdata is not None) and (self.diffFCA is not None):
            index = self.refdata[:,0]<1020
            index += self.refdata[:,0]>1200
            self.poptFCA= np.polyfit(1/self.refdata[index,0], self.diffFCA[index], 2)


    def substractFCA(self):
        if (self.diffFCA is not None) and (self.poptFCA is not None) and (self.refdata is not None):
            self.PurespectrumFCA=self.diffFCA-np.polyval(self.poptFCA,1/self.refdata[:,0])


    def CalculateOxygen(self):
        def Calalpha(x,T):
            a = -1/x*((0.09-np.exp(1.7*x))+np.sqrt((0.09-np.exp(1.7*x))**2+0.36*T**2*np.exp(1.7*x)))/0.18/T
            return a
        if self.sample_thick is not None:
            if self.ChooseSL.isChecked() and (self.diffOLine is not None):
                spectrum = self.diffOLine
            elif self.ChooseFCA.isChecked() and (self.PurespectrumFCA is not None):
                spectrum = self.PurespectrumFCA
            ind1040 = (np.abs(self.refdata[:,0]-1040)).argmin()
            ind1160 = (np.abs(self.refdata[:,0]-1160)).argmin()
            OxygenSpectrum = spectrum[ind1040:ind1160]
            OxygenWn = self.refdata[ind1040:ind1160,0]
            indpeak = OxygenSpectrum.argmax()
            Tp = 10**(-OxygenSpectrum[indpeak])
            Tb = 10** (-(OxygenSpectrum[0]+(OxygenSpectrum[-1]-OxygenSpectrum[0])/(OxygenWn[-1]-OxygenWn[0])*(OxygenWn[indpeak]-OxygenWn[0])))
            alphap = Calalpha(1e-4*self.sample_thick,Tp)
            alphab = Calalpha(1e-4*self.sample_thick,Tb)
            alphaO = alphap - alphab
            OxygenConcentration = 6.28*alphaO
            self.ax2.plot([OxygenWn[0],OxygenWn[-1]],[OxygenSpectrum[0],OxygenSpectrum[-1]],'k--')
            self.ax2.plot([OxygenWn[indpeak],OxygenWn[indpeak]],[-np.log10(Tb),OxygenSpectrum[indpeak]],'k--')
            self.ax2.plot(OxygenWn[indpeak],OxygenSpectrum[indpeak],'b*')
            self.ax2.set_title('Oxygen concentration is ' + str(round(OxygenConcentration,3)) + r' $ppma$')
            self.canvas2.draw()


    def CalculateCarbon(self):
        if self.sample_thick is not None:
            if self.ChooseSL.isChecked() and (self.diffCLine is not None):
                spectrum = self.diffCLine
            elif self.ChooseFCA.isChecked() and (self.PurespectrumFCA is not None):
                spectrum = self.PurespectrumFCA

            ind590 = (np.abs(self.refdata[:,0]-590)).argmin()
            ind620 = (np.abs(self.refdata[:,0]-620)).argmin()
            CarbonSpectrum = spectrum[ind590:ind620]
            CarbonWn = self.refdata[ind590:ind620,0]
            indpeak = CarbonSpectrum.argmax()
            Ap = CarbonSpectrum[indpeak]
            Ab = (CarbonSpectrum[0]+(CarbonSpectrum[-1]-CarbonSpectrum[0])/(CarbonWn[-1]-CarbonWn[0])*(CarbonWn[indpeak]-CarbonWn[0]))
            alphaC = 23.03*(Ap - Ab)/self.sample_thick
            CarbonConcentration = 1.64*alphaC
            self.ax3.plot([CarbonWn[0],CarbonWn[-1]],[CarbonSpectrum[0],CarbonSpectrum[-1]],'k--')
            self.ax3.plot([CarbonWn[indpeak],CarbonWn[indpeak]],[Ab,CarbonSpectrum[indpeak]],'k--')
            self.ax3.plot(CarbonWn[indpeak],CarbonSpectrum[indpeak],'b*')
            self.ax3.set_title('Carbon concentration is ' + str(round(CarbonConcentration,3)) + r' $ppma$')
            self.canvas2.draw()

    def plotdata(self):
        self.ax1.clear()
        self.ax1.set_xlabel(r'Wavenumber $[cm^{-1}]$')
        self.ax1.set_ylabel('Absorbance')
        self.ax1.set_xlim([1800,400])
        self.ax2.clear()
        self.ax2.set_title('Oxygen')
        self.ax2.set_xlim([1250,1000])
        self.ax3.clear()
        self.ax3.set_title('Carbon')
        self.ax3.set_xlabel(r'Wavenumber $[cm^{-1}]$')
        self.ax3.set_xlim([700,500])
        if self.sampledata is not None:
            self.ax1.plot(self.sampledata[:,0],self.sampledata[:,1],'r',label='Sample')
        if self.refdata is not None:
            self.ax1.plot(self.refdata[:,0],self.refdata[:,1],'b',label='Reference')
        if (self.refdatasubO is not None) and self.ChooseSL.isChecked():
            index1000 = (np.abs(self.refdata[:,0]-1000)).argmin()
            index1250 = (np.abs(self.refdata[:,0]-1250)).argmin()
            self.ax2.plot(self.refdata[index1000:index1250,0],self.refdatasubO[index1000:index1250],'b',label='Reference line subtracted')
        if (self.sampledatasubO is not None) and self.ChooseSL.isChecked():
            index1000 = (np.abs(self.sampledata[:,0]-1000)).argmin()
            index1250 = (np.abs(self.sampledata[:,0]-1250)).argmin()
            self.ax2.plot(self.sampledata[index1000:index1250,0],self.sampledatasubO[index1000:index1250],'r',label='Sample line subtracted')
        if (self.refdatasubC is not None) and self.ChooseSL.isChecked():
            index500 = (np.abs(self.refdata[:,0]-500)).argmin()
            index700 = (np.abs(self.refdata[:,0]-700)).argmin()
            self.ax3.plot(self.refdata[index500:index700,0],self.refdatasubC[index500:index700],'b',label='Reference line subtracted')
        if (self.sampledatasubC is not None) and self.ChooseSL.isChecked():
            index500 = (np.abs(self.sampledata[:,0]-500)).argmin()
            index700 = (np.abs(self.sampledata[:,0]-700)).argmin()
            self.ax3.plot(self.sampledata[index500:index700,0],self.sampledatasubC[index500:index700],'r',label='Sample line subtracted')
        if (self.diffFCA is not None) and self.ChooseFCA.isChecked():
            self.ax1.plot(self.refdata[:,0],self.diffFCA,'.',color='lime',label='Phonon substracted')
        if (self.poptFCA is not None) and self.ChooseFCA.isChecked():
            self.ax1.plot(self.refdata[:,0],np.poly1d(self.poptFCA)(1/self.refdata[:,0]),color='darkgreen',label='Fitting FCA')
        if (self.PurespectrumFCA is not None) and self.ChooseFCA.isChecked():
            index1000 = (np.abs(self.refdata[:,0]-1000)).argmin()
            index1250 = (np.abs(self.refdata[:,0]-1250)).argmin()
            index500 = (np.abs(self.refdata[:,0]-500)).argmin()
            index700 = (np.abs(self.refdata[:,0]-700)).argmin()
            self.ax1.plot(self.refdata[:,0],self.PurespectrumFCA,'k',label='FCA substracted')
            self.ax2.plot(self.refdata[index1000:index1250,0],self.PurespectrumFCA[index1000:index1250],'k',label='FCA substracted')
            self.ax3.plot(self.refdata[index500:index700,0],self.PurespectrumFCA[index500:index700],'k',label='FCA substracted')
        if (self.diffOLine is not None) and self.ChooseSL.isChecked():
            index1000 = (np.abs(self.refdata[:,0]-1000)).argmin()
            index1250 = (np.abs(self.refdata[:,0]-1250)).argmin()
            self.ax2.plot(self.refdata[index1000:index1250,0],self.diffOLine[index1000:index1250],'k',label='Phonon substracted')
        if (self.diffCLine is not None) and self.ChooseSL.isChecked():
            index500 = (np.abs(self.refdata[:,0]-500)).argmin()
            index700 = (np.abs(self.refdata[:,0]-700)).argmin()
            self.ax3.plot(self.refdata[index500:index700,0],self.diffCLine[index500:index700],'k',label='Phonon substracted')
        self.ax1.legend(loc=0)
        self.ax2.legend(loc=0)
        self.ax3.legend(loc=0)
        self.canvas1.draw()
        self.canvas2.draw()


if __name__ == '__main__':

    app = QApplication(sys.argv)
    FTIR = FTIRAnalysis()
    sys.exit(app.exec_())
