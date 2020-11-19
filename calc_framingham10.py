#!/usr/bin/python

import os
import sys
import argparse
import logging
import pandas
import numpy
import multiprocessing

__author__ = 'jmrodriguezc'

class calculate:
    '''
    Extract the correlation values
    '''
    NUM_CPUs = 10

    def __init__(self, i):
        # handle I/O files
        self.infile = i
        # extract input data ( as dataframe )
        # set index with the first column
        self.df = pandas.read_excel(self.infile, na_values=['NA'])
        # map the lowering function to all column names
        self.df.columns = map(str.lower, self.df.columns)
        # set index with the first columns
        c = str(self.df.columns[0])
        self.df = self.df.set_index(c)

    # multiprocessing methods in dataframes
    def _apply_df(self,args):
        df, func, kwargs = args
        return df.apply(func, **kwargs)
    def _apply_by_multiprocessing(self, df, func, **kwargs):
        workers = kwargs.pop('workers')
        pool = multiprocessing.Pool(processes=workers)
        result = pool.map(self._apply_df, [(d, func, kwargs)
                for d in numpy.array_split(df, workers)])
        pool.close()
        return pandas.concat(list(result))


    def _calc_f10_rcor(self, r):
        '''
        Calculate Framingham score for 10 years for each row
        '''
        # AGE integer,
        # SEXO char(1),
        # COLTOT integer,
        # HDL integer,
        # presis integer,
        # presdi integer,
        # DIABETES integer,
        # SMOKER integer,
        # VERSION_EXPERIMENTAL integer
        # 
        # Codigo_Proteo	Codigo_externo	SEQN	N_TERRITORIES	PESA Score	PESA Score_Num	Glucosa	LPA	Calcio Score	FUMADOR_RF	TABAQ	age	sexo	coltot	hdl	presis	presdi	diabetes	smoker	MED_DIABETES	MED_DISLIPEMIA	MED_HIPERTENSION	MED_OTRAS	DIABETES_RF	HIPERTENSION_RF	DISLIPEMIA_RF	Plaque thickness

        age      = r['age']
        sexo     = r['sexo']
        coltot   = r['coltot']
        hdl      = r['hdl']
        presis   = r['presis']
        presdi   = r['presdi']
        diabetes = r['diabetes']
        smoker   = r['smoker']
        # v_exp    = r['v_exp']
        # if (v_exp = 1):
        #     if (age > 74):
        #         age = 74 
        #     if (age < 30):
        #         age = 30

        if (age > 74):
            return None
        if (age < 30):
            return None

        # convert sex to integer
        if (sexo == 'H'):
            sex = 1
        else:
            sex = 0

        # coeficients for the model of Framinghan: Colesterol total
        c_160 = 0
        c200_239 = 0
        c240_279 = 0
        c280_ = 0
        if (coltot >=50 and coltot <160):
            c_160 = 1
        if (coltot >= 200 and coltot <=239):
            c200_239 = 1
        if (coltot >= 240 and coltot <=279):
            c240_279 = 1
        if (coltot >= 280 and coltot < 777):
            c280_ = 1
        # coeficients for the model of Framinghan: HDL
        h_35 = 0
        h35_44 = 0
        h45_49 = 0
        h50_59 = 0
        h60_ = 0
        if (hdl >5 and hdl <35):
            h_35 = 1
        if (hdl >=35 and hdl <45):
            h35_44 = 1
        if (hdl >=45 and hdl <50):
            h45_49 = 1
        if (hdl >=50 and hdl <60):
            h50_59 = 1
        if (hdl >=60):
            h60_ = 1           
        # coeficients for the model of Framinghan:  Blood pressure
        bp_opti = 0
        bp_norm = 0
        bp_high = 0
        bp_i = 0
        bp_ii = 0
        if (presis < 120 or presdi < 80):
            bp_opti = 1
        if ((presis >= 120 and presis <= 129) or (presdi >= 80 and presdi <=84)):
            bp_norm = 1
        if ((presis >= 130 and presis <= 139) or (presdi >= 85 and presdi <=89)):
            bp_high = 1
        if ((presis >= 140 and presis <= 159) or (presdi >= 90 and presdi <=99)):
            bp_i = 1
        if (presis >= 160 or  presdi >= 100):
            bp_ii = 1
        if (bp_opti == 1 and (bp_norm == 1 or bp_high == 1 or bp_i == 1 or bp_ii == 1)):
            bp_opti = 0
        if (bp_high == 1 and (bp_i == 1 or bp_ii == 1)):
            bp_high = 0
        if (bp_i == 1 and bp_ii == 1):
            bp_i = 0

        # calculate the Framingham score with the original parameters
        #  men
        if (sex == 1):
            l_chol =  (0.04826 * age) \
                    - (0.65945 * c_160) + (0.17692 * c200_239) + (0.50539 * c240_279) + (0.65713 * c280_) \
                    + (0.49744 * h_35) + (0.24310 * h35_44) - (0.05107 * h50_59) - (0.48660 * h60_) \
                    - (0.00226 * bp_opti) + (0.28320 *  bp_high) + (0.52168 * bp_i) + (0.61859 * bp_ii) \
                    + (0.42839 * diabetes) + (0.52337 * smoker)
        # women
        else:
            l_chol =  (0.33766 * age) - (0.00268 * (age * age)) \
                    - (0.26138 * c_160) + (0.20771 * c200_239) + (0.24385 * c240_279) + (0.53513 * c280_) \
                    + (0.84312 * h_35) + (0.37796 * h35_44) + (0.19785 * h45_49) - (0.42951 * h60_) \
                    - (0.53363 * bp_opti) - (0.06773 *  bp_high) + (0.26288 * bp_i) + (0.46573 * bp_ii) \
                    + (0.59626 * diabetes) + (0.29246 * smoker)
        
        # calculate the measure of G, which evalutes the mean values for the experiments. It is diferent ...
        # for men
        if (sex == 1):
            g_chol = 3.0975
        # for women
        if (sex == 0):
            g_chol = 9.9245
        b_chol = numpy.exp(l_chol - g_chol)
        if (sex == 1):
            framingham =  (1 - pow(0.90015, b_chol))
        if (sex == 0):
            framingham =  (1 - pow(0.96246, b_chol))

        # calculate the score but with the adapted del REGICOR which includes the GIRONA caracteristics...
        # for men
        if (sex == 1):
            g_chol = 3.489
        # for women
        if (sex == 0):
            g_chol = 10.279
        b_chol = numpy.exp(l_chol - g_chol)
        if (sex == 1):
            regicor =  (1- pow(0.951, b_chol))
        if (sex == 0):
            regicor =  (1- pow(0.978, b_chol))

        r['framingham'] = framingham
        r['regicor'] = regicor
        return r

    def calc_framingham10(self):
        '''
        Calculate Framingham score for 10 years
        '''
        df = self.df

        # calculate framingham 10 and regicor for each row
        self.df = self._apply_by_multiprocessing(df, self._calc_f10_rcor, axis=1, workers=self.NUM_CPUs)        

    
    def to_csv(self, outfile):
        '''
        Print to CSV
        '''
        if not self.df.empty:
            self.df.to_csv(outfile)
        else:
            logging.error("Empty output")

def main(args):
    ''' Main function'''

    logging.info('create calculator object')
    w = calculate(args.infile)

    logging.info('print the correlation file')
    w.calc_framingham10()

    logging.info('print the correlation file')
    w.to_csv(args.outfile)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Add the Framingham scores to input file',
        epilog='''
        Example:
            calc_framingham10.py -i data/AWHS_params_framingham.xlsx -o data/AWHS_params_framingham.sc.csv
        
        FMY: https://www.seh-lelha.org/modelos-calculo-la-probabilidad-riesgo-cardiovascular-estudio-framingham-proyecto-score/

        ''')
    parser.add_argument('-i',  '--infile', required=True, help='Excel file with Framingham parameters')
    parser.add_argument('-o',  '--outfile', required=True, help='Output file as the input file but added new scores')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info('start '+os.path.basename(__file__))
    main(args)
    logging.info('end '+os.path.basename(__file__))