# File: e (Python 2.7)

'''
Created on Mon May 26 10:55:05 2014

@author: Vlad
'''
from __init__ import um, Q_
import xlrd as xls
import numpy as np
from incert import ucreate

def GenExcelName(row, col):
    if col < 26:
        colName = chr(col + ord('A'))
    else:
        colName = chr((col / 26 - 1) + ord('A')) + chr(col % 26 + ord('A'))
    return '%s%s' % (colName, row + 1)


class InstrErr:
    
    def __init__(self, errtype = {
        'errtype': 'absolut',
        'max': 0 }, errval = 0):
        self.errtype = errtype
        self.errval = errval

    
    def getAbsolut(self):
        if self.errtype['errtype'] == 'absolut':
            return self.errval
        if self.errtype['errtype'] == 'percent':
            return self.errtype['max'] * self.errval



class OpenXls:
    """
    clasa citeste din fisierul de masuratori xls si prelucreaza date cu 
    incertitudinile existente, in functie de eroarea aleatorie (precizie)
    si eroarea absoluta(acuratetea) luata din eroare instrumentului
    
    Indexator
    ---------
    clasa contine un indexator care poate lua valori numerice si/sau literale:
    OpenXls[1,'debit_w'] -ne da debitul de la prima masuratoare 
    OpenXls[:,'debit_w'] -ne da toate debitele
    
    Observatie
    ----------
    Aceasta clasa foloseste indexarea de la 1 nu de la zero
    """
    __dtatyxls = {'nrcrt': 7,
     'debit_w': 8,
     'tin_w': 9,
     'tout_w': 10,
     'dp_w': 12,
     'm_diaf': 13,
     'dp_diaf': 14,
     'rh_a': 16,
     'pb_a': 17,
     'debit_a': 18,
     'tin_a': 19,
     'tout_a': 25,
     'dp_a': 26}
    __dtaty = {'nrcrt': 0,
     'debit_w': 1,
     'tin_w': 2,
     'tout_w': 3,
     'dp_w': 4,
     'm_diaf': 5,
     'dp_diaf': 5,
     'rh_a': 7,
     'pb_a': 8,
     'debit_a': 9,
     'tin_a': 10,
     'tout_a': 11,
     'dp_a': 12}
    __dtatyerr = {'debit_w': InstrErr(errtype={'errtype': 'percent',
                 'max': 500}, errval=0.002),
     'tin_w': InstrErr(errval=0.1),
     'tout_w': InstrErr(errval=0.1),
     'dp_w': InstrErr(errval=0.002, errtype={'errtype': 'percent',
              'max': 500}),
     'm_diaf': InstrErr(),
     'dp_diaf': InstrErr(errval=0.002, errtype={'errtype': 'percent',
                 'max': 50}),
     'rh_a': InstrErr(errval=0.001),
     'pb_a': InstrErr(errval=0.002, errtype={'errtype': 'percent',
              'max': 10}),
     'debit_a': InstrErr(errval=0.01),
     'tin_a': InstrErr(errval=0.1),
     'tout_a': InstrErr(errval=0.1),
     'dp_a': InstrErr(errval=0.1)}
    __dtaum = {'nrcrt': um.parse_expression('dimensionless'),
     'debit_w': um.parse_expression('l/min'),
     'tin_w': um.parse_expression('degC'),
     'tout_w': um.parse_expression('degC'),
     'dp_w': um.parse_expression('mbar'),
     'm_diaf': um.parse_expression('dimensionless'),
     'dp_diaf': um.parse_expression('mbar'),
     'rh_a': um.parse_expression('dimensionless'),
     'pb_a': um.parse_expression('bar'),
     'debit_a': um.parse_expression('kg/s'),
     'tin_a': um.parse_expression('degC'),
     'tout_a': um.parse_expression('degC'),
     'dp_a': um.parse_expression('mh2o')}
    
    def __skeys(self, d):
        dtatylst = sorted(d.items(), key = lambda (k, v):v)
        return map(lambda (k, v): k, dtatylst)

    
    def keys(self):
        return self._OpenXls__skeys(self._OpenXls__dtaty)

    
    def __iter__(self):
        return iter(self.udata)

    
    def __getitem__(self, pctind):
        """
        intoarce datele incarcate.
        
        Parametrii
        ----------
        pcind: int or tuple
            indexul masuratorii.
        
        Intoarce
        --------
        out: numpy.ndarray of pint.Quantity
            intoarce masuratorile cerute in functie de intrare.
            
        observatii:
        -----------
        
        ===========     ==============      ============================================      
        pcind           tip                 tipul datelor la intoarere
        ===========     ==============      ============================================
        ``:``           slice               [Q_([ufloat]),Q_([ufloat]),Q_([ufloat]),...]
        ``(x,str)``     tuple               Q_([ufloat][x])
        ===========     ==============      ============================================               
        """
        if isinstance(pctind, slice):
            return [ qu[pctind] for qu in self.udata ]
        if pctind > self.nrpcts - 1 or pctind < 0:
            Exception('index-ul este mai mare decat numarul de valori sau mai mic decat 0')
        if isinstance(pctind, tuple):
            x, y = pctind
            if isinstance(y, str):
                return self.udata[self.__dtaty[y]][x]
            if isinstance(y, int):
                return self.udata[y][x]
            if isinstance(y, slice):
                return [ qu[x] for qu in self.udata[y] ]
        if isinstance(pctind, str):
            return self.udata[self.__dtaty[pctind]]
        else:
            return [ qu[pctind] for qu in self.udata ]

    
    def open(self, fileName, sheetindex = 0, firstrow = 8, lastrow = 201):
        workbook = xls.open_workbook(fileName)
        sheet = workbook.sheet_by_index(sheetindex)
        cols = self.__dtatyxls.values()
        cols.sort() 
        self.rdata = np.array([ [ sheet.cell_value(r, c) for c in cols ] for r in range(firstrow - 1, lastrow - 1) ])
        self.nrpcts = int(self.rdata[len(self.rdata) - 1, 0])
        upcts = []
        skeys = self.__skeys(self.__dtaty)
        for i in range(1, self.nrpcts + 1):
            pct = self.rdata[self.rdata[:, 0] == i, :]
            upct = [ ucreate(pct[:, self.__dtaty[key]], self.__dtatyerr[key].getAbsolut(), self.confid) for key in skeys[1:] ]
            upct.insert(0, i)
            upcts.append(upct)

        qupcts = [ Q_(np.array(upcts)[:, self.__dtaty[key]], self.__dtaum[key]) for key in skeys[:] ]
        self.udata = qupcts
        sc = {'Lungime': sheet.cell_value(0, 38) * um.mm,
         'Latime': sheet.cell_value(0, 40) * um.mm,
         'Grosime': sheet.cell_value(0, 45) * um.mm,
         'Inalt_a': sheet.cell_value(2, 41) * um.mm,
         'Pas_a': sheet.cell_value(2, 43) * um.mm,
         'Gros_a': sheet.cell_value(3, 44) * um.mm,
         'FinTip_a': sheet.cell_value(3, 41),
         'PereteDesp': sheet.cell_value(3, 47) * um.mm,
         'Ac_w': sheet.cell_value(0, 51) * um.m ** 2,
         'At_w': sheet.cell_value(1, 51) * um.m ** 2,
         'Afr_w': sheet.cell_value(2, 51) * um.m ** 2,
         'Ac_a': sheet.cell_value(3, 51) * um.m ** 2,
         'At_a': sheet.cell_value(4, 51) * um.m ** 2,
         'Afr_a': sheet.cell_value(5, 51) * um.m ** 2,
         'Nr_a': sheet.cell_value(2, 45)}
        sc.update(Dh_a=4 * sc['Ac_a'] * sc['Grosime'] / sc['At_a'], Dh_w=4 * sc['Ac_w'] * sc['Lungime'] / sc['At_w'])
        self.__dict__.update(**sc)

    
    def __init__(self, fileName = '', sheetindex = 0, firstrow = 8, lastrow = 201, confid = 95.0):
        if not isinstance(fileName, str):
            Exception('fileName trebuie sa fie string')
        if not fileName.strip():
            self.confid = confid
            self.nrpcts = 0
            return 
        self.confid = confid
        self.open(fileName, sheetindex, firstrow, lastrow)



def testOpenXls():
    file_location = 'E:/Documents/Lucrare Doctorat/Experimental/RA 28788-0.xls'
    data = OpenXls(file_location)
    print 'running tests:.....'
    print 'data[1]\n', data[1]
    print "data[1,'debit_w']\n", data[(1, 'debit_w')]
    print 'data[1,:]\n', data[1, :]
    print "data[:,'debit_w']\n", data[:, 'debit_w'], 'count', data.nrpcts, '==', len(data[:, 'debit_w'])
    print "data['debit_w']\n", data['debit_w']
    print 'Ac_w=', str(data.Ac_w).decode('utf8').encode('ascii','ignore')
    print 'At_w=', data.At_w
    print 'Ac_a=', data.Ac_a
    print 'At_a=', data.At_a
    print 'Dh_w=', data.Dh_w
    print 'Dh_a=', data.Dh_a

if __name__ == '__main__':
    testOpenXls()
