import numpy as np
from MLG.name import name, name_good 
from MLG.Math import percentile
from MLG import  paperpath, path
from astropy.table import Table



def create_Table(Analysis_result, Analysis_result5, Analysis_result_external = None, sort_epoch =False):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	tab = None
	lines = []
	epoch = []
	rel_err = []
	lens5 = np.array([i[0].getId() for i in Analysis_result5[2]])
	source5 = np.array([i[1].getId() for i in Analysis_result5[2]])
	if Analysis_result_external is not None:
		lens_ext = np.array([i[0].getId() for i in Analysis_result_external[2]])
		source_ext = np.array([i[1].getId() for i in Analysis_result_external[2]])
		external = True
	else:
		external = False
	for i in range(len(Analysis_result[0])):
		lens = Analysis_result[2][i][0]
		s_ID1 = str(lens.getId())
		s_name = name(lens.getId(), False, True)
		Good_ = name_good(lens.getId())

		if Good_ and lens.getMag() > 5:
			s_count = '%s'%i
			source = Analysis_result[2][i][1]
			s_ID2 = str(source.getId())
			s_TCA = '%.3f'%lens.getTca()
			s_Mass = '%.2g'%lens.getMass()

			if source.getPx() == 0: fivepar = '^{*}'
			else: fivepar =''
			mm_10 = np.array([Analysis_result[0][i][x][0][6] for x in range(len(Analysis_result[0][i]))\
				if Analysis_result[0][i][x][0][4] != -999])
			mp_10 = np.array([Analysis_result[0][i][x][0][7] for x in range(len(Analysis_result[0][i]))\
				if Analysis_result[0][i][x][0][4] != -999])
			mmm_10 = percentile(mm_10)
			ppp_10 = percentile(mp_10)
			delta_M_10_p = max(ppp_10[4],-mmm_10[3])
			delta_M_10_m = min(ppp_10[3],-mmm_10[4])  
			delta_M_10 = max(ppp_10[1], -mmm_10[1])
			c=0
			while 10**(-c)>= delta_M_10 and c < 10:
					c+=1
			c+=1
			s_delta_M_10_m = str(round(delta_M_10_m-0.49999999*10**(-c),c))	
			s_delta_M_10_p = str(round(delta_M_10_p+0.49999999*10**(-c),c))
			s_delta_M_10 = str(round(delta_M_10+0.49999999*10**(-c),c))
			while(len(s_delta_M_10_m) <= c):  s_delta_M_10_m=s_delta_M_10_m+'0'
			while(s_delta_M_10_m[-c-1]) != '.': s_delta_M_10_m=s_delta_M_10_m+'0'
			while(len(s_delta_M_10_p) <= c):  s_delta_M_10_p=s_delta_M_10_p+'0'
			while(s_delta_M_10_p[-c-1]) != '.': s_delta_M_10_p=s_delta_M_10_p+'0'
			while(len(s_delta_M_10) <= c):  s_delta_M_10=s_delta_M_10+'0'
			while(s_delta_M_10[-c-1]) != '.': s_delta_M_10=s_delta_M_10+'0'
			s_delta_M_percent = '%.2g'%(delta_M_10/lens.getMass()*100)
			s_delta_M_5_m = 'NONE'
			s_delta_M_5_p = 'NONE'
			s_delta_M_5 = 'NONE'	
			delta_M_5_test = 0
			delta_M_5_p = -999
			delta_M_5_m = -999
			delta_M_5 = -999
			which5 =  np.where((lens5 == lens.getId()) & (source5 == source.getId()))[0]
			if len(which5)  == 1 : 
				k = which5[0]
				mm_5 = np.array([Analysis_result5[0][k][x][0][6] for x in range(len(Analysis_result5[0][k]))\
					if Analysis_result5[0][k][x][0][4] != -999])
				mp_5 = np.array([Analysis_result5[0][k][x][0][7] for x in range(len(Analysis_result5[0][k]))\
					if Analysis_result5[0][k][x][0][4] != -999])
				mmm_5 = percentile(mm_5)
				ppp_5 = percentile(mp_5)
				delta_M_5_p = max(ppp_5[4],-mmm_5[3])
				delta_M_5_m = min(ppp_5[3],-mmm_5[4])  
				delta_M_5 = max(ppp_5[1], -mmm_5[1])
				delta_M_5_test = delta_M_5
				c=0
				while 10**(-c)>= delta_M_5 and c < 10:
						c+=1
				c+=1
				s_delta_M_5_m = str(round(delta_M_5_m-0.49999999*10**(-c),c))	
				s_delta_M_5_p = str(round(delta_M_5_p+0.49999999*10**(-c),c))
				s_delta_M_5 = str(round(delta_M_5+0.49999999*10**(-c),c))
				while(len(s_delta_M_5_m) <= c):  s_delta_M_5_m=s_delta_M_5_m+'0'
				while(s_delta_M_5_m[-c-1]) != '.': s_delta_M_5_m=s_delta_M_5_m+'0'
				while(len(s_delta_M_5_p) <= c):  s_delta_M_5_p=s_delta_M_5_p+'0'
				while(s_delta_M_5_p[-c-1]) != '.': s_delta_M_5_p=s_delta_M_5_p+'0'
				while(len(s_delta_M_5) <= c):  s_delta_M_5=s_delta_M_5+'0'
				while(s_delta_M_5[-c-1]) != '.': s_delta_M_5=s_delta_M_5+'0'
			

			if external:
				s_delta_M_ext_m = 'NONE'
				s_delta_M_ext_p = 'NONE'
				s_delta_M_ext = 'NONE'	
				delta_M_ext_test = 0
				which_ext =  np.where((lens_ext == lens.getId()) & (source_ext == source.getId()))[0]
				if len(which_ext)  == 1 : 
					k = which_ext[0]
					mm_ext = np.array([Analysis_result_external[0][k][x][0][6]\
						for x in range(len(Analysis_result_external[0][k])) \
						if Analysis_result_external[0][k][x][0][4] != -999])
					mp_ext = np.array([Analysis_result_external[0][k][x][0][7] \
						for x in range(len(Analysis_result_external[0][k])) \
						if Analysis_result_external[0][k][x][0][4] != -999])
					mmm_ext = percentile(mm_ext)
					ppp_ext = percentile(mp_ext)
					delta_M_ext_p = max(ppp_ext[4],-mmm_ext[3])
					delta_M_ext_m = min(ppp_ext[3],-mmm_ext[4])  
					delta_M_ext = max(ppp_ext[1], -mmm_ext[1])
					delta_M_ext_test = delta_M_ext
					c=0
					while 10**(-c)>= delta_M_ext and c < 10:
							c+=1
					c+=1
					s_delta_M_ext_m = str(round(delta_M_ext_m-0.49999999*10**(-c),c))	
					s_delta_M_ext_p = str(round(delta_M_ext_p+0.49999999*10**(-c),c))
					s_delta_M_ext = str(round(delta_M_ext+0.49999999*10**(-c),c))
					while(len(s_delta_M_ext_m) <= c):  s_delta_M_ext_m=s_delta_M_ext_m+'0'
					while(s_delta_M_ext_m[-c-1]) != '.': s_delta_M_ext_m=s_delta_M_ext_m+'0'
					while(len(s_delta_M_ext_p) <= c):  s_delta_M_ext_p=s_delta_M_ext_p+'0'
					while(s_delta_M_ext_p[-c-1]) != '.': s_delta_M_ext_p=s_delta_M_ext_p+'0'
					while(len(s_delta_M_ext) <= c):  s_delta_M_ext=s_delta_M_ext+'0'
					while(s_delta_M_ext[-c-1]) != '.': s_delta_M_ext=s_delta_M_ext+'0'
					s_delta_M_percent_ext = '%.2g'%(delta_M_ext/lens.getMass()*100)

				#---------------------------
				#create Astropy.table
				if tab is None:
					tab  = Table(names = ['name', 'lens_id', 'source_id', 'fivepar', 'TCA', 'Mass',\
						'deltaM_10','sigma_deltaM_10_+', 'sigma_deltaM_10_-',\
						'deltaM_5' , 'sigma_deltaM_5_+','sigma_deltaM_5_-',\
						'deltaM_ext' , 'sigma_deltaM_ext_+','sigma_deltaM_ext_-'],\
						dtype = [object, np.int64, np.int64, np.bool_,np.float64, np.float64, np.float64, np.float64,\
						np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64]) 
				tab.add_row([s_name,lens.getId(),source.getId(),source.getPx() != 0,lens.getTca(),lens.getMass(),\
					delta_M_10, delta_M_10_p, delta_M_10_m,\
					delta_M_5, delta_M_5_p, delta_M_5_m,\
					delta_M_ext, delta_M_ext_p, delta_M_ext_m])			
				#---------------------------

				#---------------------------
				#element list for sorting 
				epoch.append(lens.getTca())
				rel_err.append(delta_M_10/lens.getMass())
				#---------------------------

				#---------------------------
				#create row in latex table 
				# \\(int\\) & \\(Name\\) & \\(ID1\\) & \\(ID2\\) & \\(T_CA\\) & \\(M_in\\)
				# & \\(deltaM_plus^{sigma+}_{sigma-}\\) & \\(deltaM_minus^{sigma+}_{sigma-}  \\) \\\\
				line = '%s & \\(%s\\) & \\(%s%s\\) & \\(%s\\) & \\(%s\\) & \\(\\pm%s^{+%s}_{%s}\\)  & \\(%s%s\\)'\
						%(s_name,s_ID1, s_ID2,fivepar, s_TCA, s_Mass,\
							s_delta_M_10, s_delta_M_10_p, s_delta_M_10_m, s_delta_M_percent,'\\%')
				if s_delta_M_5_m == 'NONE' or  delta_M_5_test > lens.getMass():
					line = line + ' & '
				else:
					line = line + ' & \\(\\pm%s^{+%s}_{%s}\\)'% (s_delta_M_5, s_delta_M_5_p, s_delta_M_5_m)
				if s_delta_M_ext_m == 'NONE' or  delta_M_ext_test > lens.getMass():
					line = line + ' & '
				else:
					
					line = line + ' & \\(%s%s\\)'%(s_delta_M_percent_ext,'\\%')
				lines.append(line) 

				#---------------------------
				
			else:
				#---------------------------
				#create Astropy.table
				if tab is None:
					tab  = Table(names = ['name', 'lens_id', 'source_id', 'fivepar', 'TCA', 'Mass',\
							'deltaM_10', 'sigma_deltaM_10_+', 'sigma_deltaM_10_-',\
							'deltaM_5' , 'sigma_deltaM_5_+','sigma_deltaM_5_-'],\
							dtype = [object, np.int64, np.int64, np.bool_,np.float64, np.float64, np.float64,\
							np.float64, np.float64, np.float64, np.float64, np.float64]) 
				tab.add_row([s_name,lens.getId(),source.getId(),source.getPx() != 0,lens.getTca(),lens.getMass(),\
					delta_M_10, delta_M_10_p, delta_M_10_m, delta_M_5, delta_M_5_p, delta_M_5_m])			
				#---------------------------
				#element in list for sorting 
				epoch.append(lens.getTca())
				rel_err.append(delta_M_10/lens.getMass())
				#---------------------------

				#---------------------------
				#create row in latex table 
				#'\\(int\\) & \\(Name\\) & \\(ID1\\) & \\(ID2\\) & \\(T_CA\\) & \\(M_in\\) 
				#& \\(deltaM_plus^{sigma+}_{sigma-}\\) & \\(deltaM_minus^{sigma+}_{sigma-}  \\) \\\\'
				line = '%s & \\(%s\\) & \\(%s%s\\) & \\(%s\\) & \\(%s\\) & \\(\\pm%s^{+%s}_{%s}\\)  & \\(%s%s\\)'\
						%(s_name,s_ID1, s_ID2,fivepar, s_TCA, s_Mass, s_delta_M_10, s_delta_M_10_p, s_delta_M_10_m, s_delta_M_percent,'\\%')
				if s_delta_M_5_m == 'NONE' or  delta_M_5_test > lens.getMass():
					line = line + ' & '
				else:
					line = line + ' & \\(\\pm%s^{+%s}_{%s}\\)'% (s_delta_M_5, s_delta_M_5_p, s_delta_M_5_m)
				line = line + ' \\\\'
				lines.append(line) 
	
	#---------------------------
	#sorting list		
	if sort_epoch:	
		lines = [x for _,x in sorted(zip(epoch,lines))]
		rel_err = [x for _,x in sorted(zip(epoch,rel_err))]
		epoch = [x for _,x in sorted(zip(epoch,epoch))]
	else:
		lines = [x for _,x in sorted(zip(rel_err,lines))]
		epoch = [x for _,x in sorted(zip(rel_err,epoch))]
		rel_err = [x for _,x in sorted(zip(rel_err,rel_err))]
	#---------------------------

	#---------------------------
	#save Table
	if external: tablename = 'result_single_ext'
	else: tablename = 'result_single'
	print('write table: ' + tablename+ '.vot')
	print('write tex_table: ' + tablename+ '.txt')
	tab.write(path + 'Data/' + tablename+ '.vot', format = 'votable',overwrite=True)
	#---------------------------
	
	#---------------------------
	#writing Latex Table
	f = open(paperpath+tablename+'.tex','w')
	
	#Setup a Table 
	f.write('\\renewcommand*{\\arraystretch}{1.4}'+'\n')
	f.write('\\begin{table*}[]' +'\n')
	f.write('\\input{result_single_caption1}\n') # include caption
	f.write('\\label{table:single}' +'\n')
	f.write('\\tiny')
	
	#Format and head Row
	f.write('\\begin{tabular}{l|rrrrr|rr|r|}' +'\n')
	f.write('\\# & Name-Lens & \\(DR2\\_ID\\)-Lens & \\(DR2\\_ID\\)-Source & \\(T_{CA}\\) & \\(M_{in}\\)')
	f.write(' & \\(\\sigma\\,M_{10}\\) & \\(\\sigma\\,M_{10}/M_{in}\\) & \\(\\sigma\\,M_{5}  \\)')
	f.write(' \\\\\n')
	f.write('  & & & & \\(\\mathrm{Jyear}\\) & \\(\\mathrm{M_{\\odot}}\\)')
	f.write(' & \\(\\mathrm{M_{\\odot}}\\) & & \\(\\mathrm{M_{\\odot}}\\)')
	#if external: f.write(' & \\(\\mathrm{M_{\\odot}}\\) \\\\\n')
	f.write(' \\\\\n')
	f.write('\\hline' +'\n')
	
	#
	c3 = 1 #row counter combined
	c4 = 0 #row counter individual tables


	if sort_epoch: #sorted by epoch
		#write rows
		for l in range(len(lines)):
			if (epoch[l] <= 2019.5) & (rel_err[l] <=0.5): 
				if c4%10> 4: f.write('\\rowcolor{lightGray}\n')	
				ll = lines[l]
				if external:
					ll = ll.split('&')
					ll.pop(-1)
					ll = '&'.join(ll)
				f.write('\\(%d\\) & '%c3 + ll+'\\\\\n')
				c3 +=1
				c4+=1
		#Table end		
		f.write('\\hline' +'\n')
		f.write('\\end{tabular}' + '\n' + '\\end{table*}')

		#Setup second Table
		f.write('\n'+'\n'+'\n'+'\n'+'\n'+'\\begin{table*}[]' +'\n')
		f.write('\\input{result_single_caption2}\n')
		f.write('\\label{table:single2}' +'\n')
		f.write('\\tiny')

		#Format and head Row
		f.write('\\begin{tabular}{l|rrrrr|rr|r|}' +'\n')
		f.write('\\# & Name-Lens & \\(DR2\\_ID\\)-Lens & \\(DR2\\_ID\\)-Source & \\(T_{CA}\\) & \\(M_{in}\\)')
		f.write(' & \\(\\sigma\\,M_{10}\\) & \\(\\sigma\\,M_{10}/M_{in}\\)')
		if external: f.write(' & \\(\\sigma\\,M_{obs}/M_{in}\\)\\\\\n')
		else: f.write(' & \\(\\sigma\\,M_{5}\\) \\\\\n')
		f.write('  & & & & \\(\\mathrm{Jyear}\\) & \\(\\mathrm{M_{\\odot}}\\)')
		if external: f.write(' & \\(\\mathrm{M_{\\odot}}\\) & & \\\\\n')
		else: f.write(' & \\(\\mathrm{M_{\\odot}}\\) & & \\(\\mathrm{M_{\\odot}}\\)\\\\\n')
		f.write('\\hline' +'\n')

		
		c4 = 0 #Row counter individual tables
		
		#Write rows
		for l in range(len(lines)):
			if (epoch[l] > 2019.5) & (rel_err[l] <=0.5): 
				if c4%10> 4: f.write('\\rowcolor{lightGray}\n')	
				ll=lines[l]
				if external:
					ll = ll.split('&')
					ll.pop(-2)
					ll = '&'.join(ll)
				f.write('\\(%d\\) & '%c3 + ll+'\\\\\n')
				c3 +=1
				c4+=1

		#Table end		
		f.write('\\hline' +'\n')
		f.write('\\end{tabular}' + '\n' + '\\end{table*}')	
		f.close()
	

	else:
		i = 0
		c2 = [15,30,50,100,101,200]
		for line in lines: 
			if c2[i]<= 50:
				if c4%10> 4: f.write('\\rowcolor{lightGray}\n')	
				f.write('\\(%d\\) & '%c3 + line+'\\\\\n')
				c3 +=1
				c4+=1
			if c < len(rel_err):
				if rel_err[c]*100 > c2[i]:
					while rel_err[c]*100 > c2[i]:
						i+=1  
					if c2[i] == 50: 
						f.write('\\hline' +'\n')
						f.write('\\end{tabular}' + '\n' + '\\end{table*}')
						f.write('\n'+'\n'+'\n'+'\n'+'\n'+'\\begin{table*}[]' +'\n')
						f.write('\\input{result_single_caption2}\n')
						f.write('\\label{table:single2}' +'\n')
						f.write('\\tiny')
						if external: f.write('\\begin{tabular}{l|rrrrr|rr|r|r|}' +'\n')
						else: f.write('\\begin{tabular}{l|rrrrr|rr|r|}' +'\n')
						f.write('\\# & Name-Lens & \\(DR2\\_ID\\)-Lens & \\(DR2\\_ID\\)-Source')
						f.write(' & \\(T_{CA}\\) & \\(M_{in}\\)')
						f.write(' & \\(\\sigma\\,M_{10}\\) & \\(\\sigma\\,M_{10}/M_{in}\\) & \\(\\sigma\\,M_{5}  \\)')
						if external: f.write(' & \\(\\sigma\\,M_{obs}  \\) \\\\\n')
						else: f.write(' \\\\\n')
						f.write('  & & & & \\(\\mathrm{Jyear}\\) & \\(\\mathrm{M_{\\odot}}\\) \\\\\n')
						f.write(' & \\(\\mathrm{M_{\\odot}}\\) & & \\(\\mathrm{M_{\\odot}}\\) \\\\\n')
						f.write('\\hline' +'\n')
						c4 = 0 

			if c2[i] >= 100: break
			c+=1
		f.write('\\hline' +'\n')
		f.write('\\end{tabular}' + '\n' + '\\end{table*}')	
		f.close()

def create_Table_multi(Analysis_multi_result):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	Dat_single, Dat_all, Dat_fps, Dat_sig = Analysis_multi_result
	
	q1 = []
	q2 = []



	for i in range(len(Dat_all[1])):
		lens = Dat_all[2][i][0]
		s_Mass = '%.2g'%lens.getMass()
		lens_id = Dat_all[1][i]
		s_lens_id = str(lens_id)
		s_name = name(lens.getId(), False, True)
		all_bool = True
		fps_bool = False
		sig_bool = False
		for j in range(len(Dat_fps[1])):
			if lens_id == Dat_fps[1][j]: 
				fps_bool = True
				break

		for k in range(len(Dat_sig[1])):
			if lens_id == Dat_sig[1][k]: 
				sig_bool = True
				break

		Vs_delta_M_m = ['','','']
		Vs_delta_M_p = ['','','']
		Vs_delta_M = ['','','']
		rel_err = ['','','']
		for uu in range(3):
			if [all_bool, fps_bool, sig_bool][uu]:
				nn = [i,j,k][uu]
				Analysis_result = [Dat_all, Dat_fps, Dat_sig][uu]

				mm = np.array([Analysis_result[0][nn][x][0][6] for x in range(len(Analysis_result[0][nn])) if Analysis_result[0][nn][x][0][4] != -999])
				mp = np.array([Analysis_result[0][nn][x][0][7] for x in range(len(Analysis_result[0][nn])) if Analysis_result[0][nn][x][0][4] != -999])
				mmm = percentile(mm)
				ppp = percentile(mp)
				delta_M_p = max(ppp[4],-mmm[3])
				delta_M_m = min(ppp[3],-mmm[4])  
				delta_M = max(ppp[1], -mmm[1])
				c=0
				while 10**(-c)>= delta_M and c < 10:
						c+=1
				c+=1
				s_delta_M_m = str(round(delta_M_m-0.49999999*10**(-c),c))	
				s_delta_M_p = str(round(delta_M_p+0.49999999*10**(-c),c))
				s_delta_M = str(round(delta_M+0.49999999*10**(-c),c))
				while(len(s_delta_M_m) <= c):  s_delta_M_m=s_delta_M_m+'0'
				while(s_delta_M_m[-c-1]) != '.': s_delta_M_m=s_delta_M_m+'0'
				while(len(s_delta_M_p) <= c):  s_delta_M_p=s_delta_M_p+'0'
				while(s_delta_M_p[-c-1]) != '.': s_delta_M_p=s_delta_M_p+'0'
				while(len(s_delta_M) <= c):  s_delta_M=s_delta_M+'0'
				while(s_delta_M[-c-1]) != '.': s_delta_M=s_delta_M+'0'
				Vs_delta_M_m[uu] = s_delta_M_m
				Vs_delta_M_p[uu] = s_delta_M_p
				Vs_delta_M[uu] = s_delta_M
			if uu ==0:
				rel_err1 = (delta_M/lens.getMass())
			rel_err[uu] = '%.2g'%(delta_M/lens.getMass()*100)

	# sortiere Source ID into 4 
		ID_all = [x.getId() for x in Dat_all[2][i][1:]]
		ID_fps = [x.getId() for x in Dat_fps[2][j][1:]]
		ID_sig = [x.getId() for x in Dat_sig[2][k][1:]]
		tps = []
		fps = []
		sig = []
		for ID in  ID_all:
			if ID in ID_sig: sig.append(ID)
			elif ID in ID_fps: fps.append(ID)
			else: tps.append(ID)

		q1.append([rel_err1,rel_err,s_Mass, s_lens_id,s_name, tps,fps ,sig, Vs_delta_M_m,Vs_delta_M_p,Vs_delta_M])
	q1.sort(key=lambda x: x[0])
		
	f = open(paperpath+'result_multi.tex','w')
	f.write('\\renewcommand*{\\arraystretch}{1.4}' + '\n')
	f.write('\\begin{sidewaystable*}[]' + '\n')
	f.write('\\input{result_multi_caption1}\n')
	f.write('\\label{table:multi}'+ '\n')
	f.write('\\tiny\\begin{tabular}{l|r|rrrr|rr|rr|rr|}' + '\n')
	f.write('\\# & Name-Lens & ' + '\n')
	f.write('\\(DR2\\_ID\\)-Source & \\(DR2\\_ID\\)-Source' )
	f.write(' & \\(DR2\\_ID\\)-Source' + '\n')
	f.write(' & \\(M_{\mathrm{in}}\\) & \\(\\sigma\\,M_{\mathrm{all}}\\)')
	f.write(' & \\(\\sigma\\,M_{\mathrm{all}}/M_{\mathrm{in}}\\) & \\(\\sigma\\,M_{\mathrm{5-par}}  \\)')
	f.write(' & \\(\\sigma\\,M_{\mathrm{5-par.}}/M_{in}\\)  & \\(\\sigma\\,M_{\mathrm{sig.}}  \\)')
	f.write(' & \\(\\sigma\\,M_{\mathrm{sig.}}/M_{\mathrm{in}}\\)  \\\\' + '\n')
	f.write(' & \\(DR2\\_ID\\)-Lens & (2-parameter) & (5-parameter)& (sigma)')
	f.write(' &  \\(\\Msun{}\\) & \\(\\Msun{}\\) & \\(\\%\\)& \\(\\Msun{}\\) & \\(\\%\\)& \\(\\Msun{}\\) & \\(\\%\\)')
	f.write(' \\\\\n')
	f.write('\\hline' +'\n')

	for n  in range(len(q1)):

		rel_err1,rel_err, s_Mass, s_lens_id,s_name,  tps,fps,sig,  Vs_delta_M_m,Vs_delta_M_p,Vs_delta_M =q1[n]
		if rel_err1 < 0.5:
			ii = 0
			if s_name == '' :ll = ['\\(%s\\)'% s_lens_id,]
			else: ll = [s_name, '\\(%s\\)'% s_lens_id]
			for line_index in range(max([len(tps), len(fps), len(sig),len(ll)])):
				if line_index >= len(tps): s_tps = ''
				else: s_tps = str(tps[line_index])
				if line_index >= len(fps): s_fps = ''
				else: s_fps = str(fps[line_index])
				if line_index >= len(sig): s_sig = ''
				else: s_sig = str(sig[line_index])
				if line_index >= len(ll): s_ln = ''
				else: s_ln= str(ll[line_index])

				if line_index == 0: 
					line = '%s &  %s  & \\(%s\\) & \\(%s\\) & \\(%s\\)  & \\(%s\\)  &'\
						%(str(n+1), s_ln, s_tps,s_fps, s_sig, s_Mass)
					if len(tps) == 0: line = line + ' & &' 
					else: line = line+ '\\(\\pm%s^{+%s}_{%s}\\) & \\(%s%s\\) &' \
						% (Vs_delta_M[0], Vs_delta_M_p[0], Vs_delta_M_m[0], rel_err[0], '\\%')
					if len(fps) == 0: line = line + ' & &'
					else:	
						line = line+ '\\(\\pm%s^{+%s}_{%s}\\) &\\(%s%s\\) &' \
						% (Vs_delta_M[1], Vs_delta_M_p[1], Vs_delta_M_m[1], rel_err[1],'\\%')
					if len(sig) == 0: line = line + ' \\\\'
					else:	
						line = line+ '\\(\\pm%s^{+%s}_{%s}\\) & \\(%s%s\\)  \\\\'\
							% (Vs_delta_M[2], Vs_delta_M_p[2], Vs_delta_M_m[2], rel_err[2],'\\%')
					 
				else: 
					line = ' &%s& \\(%s\\) & \\(%s\\) & \\(%s\\)  & & & & & & & \\\\'\
						%  (s_ln, s_tps,s_fps, s_sig)

				if n%2 : f.write('\\rowcolor{lightGray}\n')	
				f.write(line)
				f.write('\n')
			f.write('\\hline'+ '\n')


	f.write('\\end{tabular}' +'\n' +'\\end{sidewaystable*}')
	f.close()