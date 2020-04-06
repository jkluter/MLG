import matplotlib.pyplot as plt
import numpy  as np	
from MLG import imagepath, paperpath, path
from imageio import imread
import matplotlib.cbook as cbook
from MLG.utils import color_own 
from matplotlib import rc

__all__ = ['dark','Image_Window','Image_precision','Image_Illustration','Image_Illustration_Multi','Image_compare_micro','Image_astroshift', 'create_all_Image']
def dark(onof = 0):
	if onof is 'on': plt.style.use('dark_background')
	elif onof is 'off': plt.style.use('default')
	elif onof == True: plt.style.use('dark_background')
	else: plt.style.use('default')

def Image_Window(string = 'resolve_Window', pres = False):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
	else: black = color_own([0.,0.,0.,1])

	c1 = color_own([0,1,1,1])
	c2 = color_own([1,0,0,1])
	c3 = color_own([1,1,0.2,1])
	c4 = color_own([0.4,0.4,0.4,1])
	c_star = [c1,c2,c3,c4]
	c_grid1 = color_own([0,int(dark),1,1])
	c_grid2 = color_own([0.5,1,0,1])
	star= np.array([[0, 0],[0.1,0.9],[-1,-1.1],[-0.5,0.1]])
	fig = plt.figure(figsize = [12,12])
	x_width = 0.059
	y_width = 0.177

	#------------------------------------------------------------
	# axis 	
	plt.xticks( fontsize = 25)
	plt.yticks( fontsize = 25)
	plt.grid(True)
	plt.axis('equal')
	plt.axis([-1.5,1.5,-1.4,1.6])
	plt.ylabel('Across-scan direction (AC) [arcsec]', fontsize = 30)
	plt.xlabel('Along-scan direction (AL) [arcsec]', fontsize = 30)

	#------------------------------------------------------------
	#------------------------------------------------------------
	#Grid Major Star

	for i in range(-6,7):
		plt.plot([-6*x_width,6*x_width], [i*y_width,i*y_width], c = c_grid1,linewidth = 3)
		plt.plot([i*x_width,i*x_width], [-6*y_width,6*y_width], c = c_grid1, linewidth = 3)
	plt.text(0,1.4,"Along-scan direction\n $12\,\mathrm{pix} \\times 0.059 \mathrm{''/pix} = 0.708\mathrm{''}$",fontsize = 25, verticalalignment = 'center',  horizontalalignment = 'center', rotation = 0)	
	plt.text(0.7,0,"Across-scan direction\n $12\,\mathrm{pix} \\times 0.177 \mathrm{''/pix} = 2.124\mathrm{''}$",fontsize = 25, verticalalignment = 'center',  horizontalalignment = 'center', rotation = 90)	

	plt.arrow(0,6*y_width+2*x_width, -6*x_width+0.02,0,color= black, head_width=0.1,\
				overhang = 0.5, length_includes_head=True ,zorder = 10, linewidth = 3)
	plt.arrow(0,6*y_width+2*x_width, 6*x_width-0.02,0,color= black, head_width=0.1,\
				overhang = 0.5, length_includes_head=True ,zorder = 10, linewidth = 3)
	plt.arrow(8*x_width,0,0, -6*y_width+0.02,color= black, head_width=0.1,\
				overhang = 0.5, length_includes_head=True ,zorder = 10, linewidth = 3)
	plt.arrow(8*x_width,0,0, 6*y_width-0.02,color= black, head_width=0.1,\
				overhang = 0.5, length_includes_head=True ,zorder = 10, linewidth = 3)
	plt.scatter(star[:1,0], star[:1,1], marker=(5, 1),c = c_star[:1], s = [3000], zorder = 1000)
	if pres: fig.savefig(imagepath + string + '_1.png', format = 'png')	

	#------------------------------------------------------------
	
	#------------------------------------------------------------
	#Grid Minor Star
	plt.scatter(star[1:3,0], star[1:3,1], marker=(5, 1),c = c_star[1:3], s = [2000,2000], zorder = 1000)
	if pres: fig.savefig(imagepath + string + '_2.png', format = 'png')	
	for i in range(-5,8):
		plt.plot([-15*x_width,-6*x_width], [i*y_width,i*y_width], c = c_grid2,linewidth = 3, zorder = -1)
	for i in range(-15,-5):
		plt.plot([i*x_width,i*x_width], [-5*y_width,7*y_width], c = c_grid2, linewidth = 3, zorder = -1)
	plt.scatter(star[3:,0], star[3:,1], marker=(5, 1),c = c_star[3:], s = [2000], zorder = 1000)

	if pres: fig.savefig(imagepath + string + '_3.png', format = 'png')	

	#------------------------------------------------------------





	
	fig.savefig(imagepath + string + '.png', format = 'png')	
	print('Create Image: '+ imagepath+ string + '.png')
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	plt.close(fig)

def Image_precision(string = 'Sig_vs_Gmag', Gaia_precision = path+'InputTable/resolution_Gaia.png', pres = False):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
		color1 = color_own([0.85,0,0,1])
		color2 = color_own([0,0,1,1])
		color3 = color_own([0,1,1,1])
		color4 = color_own([0.5,1,0,1])
		color5 = color_own([1,1,0,1])
	else:
		black = color_own([0.,0.,0.,1])
		color1 = color_own([0.85,0,0,1])
		color2 = color_own([0,0,1,1])
		color3 = color_own([0,1,1,1])
		color4 = color_own([0,1,0,1])
		color5 = color_own([1,1,0,1])
	fig = plt.figure(figsize = [12,10])
	Gmag = np.arange(4,22,0.01)


	datafile = cbook.get_sample_data(Gaia_precision)
	img = imread(datafile)
	z = 10 ** (0.4 * (np.maximum(Gmag, 14) - 15)) #(14-np.minimum(Gmag, 14))
	z2 = 10 ** (0.4 * (np.maximum(Gmag, 12) - 15))
	sig_pi = (-1.631 + 680.766 * z2 + 32.732 * z2**2)**0.5/1000
	sig_fov2 =(-1.631 + 680.766 * z + 32.732 * z**2)**0.5/1000 *7.75 +0.1
	sig_fov3 =  sig_fov2 / np.sqrt(9) 
	plt.plot([0,1],[-5,-5], c = color1, linewidth = 3, label = 'formal precision from Gaia DR2 (per CCD)' )
	plt.plot([0,1],[-5,-5], c = color2, linewidth = 3, label = 'actual precision from Gaia DR2 (per CCD)' )
	plt.yticks([np.log10(i) for i in [20,10, 5,2,1, 0.5,0.2,0.1, 0.05,0.02, 0.01]],[20,10, 5,2,1, 0.5,0.2,0.1, 0.05,0.02,0.01],  fontsize = 25)
	plt.xticks( fontsize = 25)
	plt.ylabel('Standard deviation of AL field angle [mas]', fontsize = 30)
	plt.xlabel('G magnitude', fontsize = 30)	

	plt.imshow(img, zorder=0, extent=[5, 21.04, np.log10(0.0195),np.log10(10)])
	plt.axis('auto')
	plt.xlim([4,22]) 
	plt.ylim([np.log10(0.005),np.log10(40)]) 
	if pres:
		plt.legend(loc = 'upper left',fontsize = 20)
		fig.savefig(imagepath + string + '_1.png', format = 'png')	
	plt.plot(Gmag,np.log10(sig_pi), '--',c = color3, dashes =(5,5), linewidth = 3,  label= 'predicted end-of-mission parallax error')
	if pres:
		plt.legend(loc = 'upper left',fontsize = 20)
		fig.savefig(imagepath + string + '_2.png', format = 'png')	
	plt.plot(Gmag,np.log10(sig_fov2), ':' , c = color4, linewidth = 5, label= 'used Standard deviation (per CCD)' )
	if pres:
		plt.legend(loc = 'upper left',fontsize = 20)
		fig.savefig(imagepath + string + '_3.png', format = 'png')	

	plt.plot(Gmag,np.log10(sig_fov3) ,c = color5,linewidth = 7,  label= 'used Standard deviation for 9 CCD observations' )
	plt.plot([5, 21.04, 21.04,5,5], [np.log10(0.0195),np.log10(0.0195),np.log10(10),np.log10(10),np.log10(0.0195)], linewidth = 2, color = [0.5,0.5,0.5,1], zorder = 0.1)
	plt.axis('auto')
	plt.xlim([4,22]) 
	plt.ylim([np.log10(0.005),np.log10(40)]) 
	plt.legend(loc = 'upper left',fontsize = 20)

	fig.savefig(imagepath + string + '.png', format = 'png')	
	print('Create Image: '+ imagepath+ string + '.png')
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	plt.close(fig)

def Image_Illustration(string = 'Illustration'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''

	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
		color1 = color_own([0,1,1,1])
		color2 = color_own([1,0.5,0,1])
		color3 = color_own([0.5,1,0,1])
		color4 = color_own([1,0,1,1])
		color5 = color_own([0,1,1,1])
		color6 = color_own([0,1,1,1])

	else: 
		black = color_own([0.,0.,0.,1])
		color1 = color_own([0,0,1,1])
		color2 = color_own([1,0.5,0,1])
		color3 = color_own([0,1,0,1])
		color4 = color_own([1,0,1,1])
		color5 = color_own([0,1,1,1])
		color6 = color_own([0,1,1,1])

	t = np.array([12, 35, 41, 61, 73, 89])
	scandir = np.array([0.1, 0.7, 0.4, 0.8 , 0.2, 0.1])*np.pi	

	x1 = np.linspace(1,13,100) + 0.3 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y1 = np.linspace(1,3,100) + 0.3* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	

	x2 = np.linspace(3,7,100)# +  0.03 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y2 = np.linspace(7,4.5,100)# + 0.03* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	

	d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
	TE = 1.5
	X2 = x2 + (x2-x1) * TE/(d**2 +2*TE)
	Y2 = y2 + (y2-y1) * TE/(d**2 +2*TE)	

	dX2 = x1-X2
	dY2 = y1-Y2
	dx2 = x1-x2
	dy2 = y1-y2	

	fig = plt.figure(figsize= (12,8))
	ax = plt.subplot(111)
	ax.axis('equal')
	ax.axis('off')	
	
	for i in range(len(t)): 
		xm1 =np.array([-1,1]) * np.cos(scandir[i]) + x1[t[i]]
		ym1 =np.array([-1,1]) * np.sin(scandir[i]) + y1[t[i]]
		xm2 =np.array([-1,1]) * np.cos(scandir[i]) + X2[t[i]]
		ym2 =np.array([-1,1]) * np.sin(scandir[i]) + Y2[t[i]]	

		dsc = ((dx2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dy2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]
		dSC = ((dX2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dY2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]	

		ttX2 = np.array([0,-dSC[1]/2,dSC[1]/2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],X2[t[i]],X2[t[i]]])
		ttY2 = np.array([0,-dSC[1]/2,dSC[1]/2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],Y2[t[i]],Y2[t[i]]])
		ttx2 = np.array([0,-dsc[1]/2-0.2,dsc[1]/2-0.2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],x2[t[i]],x2[t[i]]])
		tty2 = np.array([0,-dsc[1]/2-0.2,dsc[1]/2-0.2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],y2[t[i]],y2[t[i]]])	

		if i % 2  == 0:
			plt.arrow(ttx2[2],tty2[2], 0.0001*(ttx2[2]-ttx2[1]),0.0001*(tty2[2]-tty2[1]),color= color1, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttx2[1],tty2[1], 0.0001*(ttx2[1]-ttx2[2]),0.0001*(tty2[1]-tty2[2]),color= color1, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[2],ttY2[2], 0.0001*(ttX2[2]-ttX2[1]),0.0001*(ttY2[2]-ttY2[1]),color= color2, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[1],ttY2[1], 0.0001*(ttX2[1]-ttX2[2]),0.0001*(ttY2[1]-ttY2[2]),color= color2, head_width=0.2,\
				overhang = 0.5, length_includes_head=True, zorder = 10)
			plt.plot(ttx2[0:2],tty2[0:2],color = black, linestyle= ':')
			plt.plot(ttX2[0:2],ttY2[0:2],color = black, linestyle= ':')
			plt.plot(ttx2[1:3],tty2[1:3],color = color1,linewidth = 3 ,  linestyle= '--',dashes=(10, 10))
			plt.plot(ttX2[1:3],ttY2[1:3],color = color2, linewidth = 3,linestyle= '-')
			plt.plot(ttx2[2:],tty2[2:],color = black, linestyle= ':')
			plt.plot(ttX2[2:],ttY2[2:],color = black, linestyle= ':')
		if i% 2 == 0: 
			plt.plot(xm2,ym2, color = black, linewidth = 3,zorder = 1)
			plt.plot(xm1,ym1, color = black, linewidth = 3,zorder = 1)
		else:  
			plt.plot(xm2,ym2, color = 'grey', linewidth = 2, zorder = -1)
			plt.plot(xm1,ym1, color = 'grey', linewidth = 2, zorder = -1)
		#if i ==0 :
	plt.plot(x1,y1, color = color3, linewidth = 3)
	plt.plot(x2,y2, color = color1, linestyle= '--',dashes=(10, 5), linewidth = 3, zorder = -1)
	plt.plot(X2,Y2, color = color2, linewidth = 3)	

	plt.xlim([-0.5,14])	

	xr = 12
	yr = 7
	plt.text(xr-0.8,0,'RA $\cdot$ cos(Dec)',verticalalignment = 'center',fontsize = 25)
	plt.text(0,yr + 0.25,'Dec',fontsize = 25, horizontalalignment = 'center', rotation = 90)	

	plt.arrow(-0.025,0,xr-1,0,width = 0.05,overhang = 0.5,head_width = 0.5, head_length = 0.5,color= black, zorder = 100,length_includes_head=True)
	plt.arrow(0,-0.025,0,yr-0.5,width = 0.05,overhang = 0.5,head_width = 0.5,head_length = 0.5,color= black, zorder = 100,length_includes_head=True)	
	
	plt.text(2,1.5,'Lens',color = color3, fontsize = 25, horizontalalignment = 'center', rotation = 0, weight = 'bold')	
	plt.text(4,7.5,'Star 1',color = color2,fontsize = 25, horizontalalignment = 'center', rotation = 0,weight = 'bold')	

	fig.savefig(imagepath + string + '.png', format = 'png')	
	print('Create Image: '+ imagepath+ string + '.png')
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	plt.close(fig)


def Image_Illustration2 (string = 'Illustration'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''

	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.,0.,0.,1])
		grey = color_own([.5,.5,0.5,1])
		cyan = color_own([0,1,1,1])
		blue = color_own([0,0,1,1])
		lime = color_own([0.6,1.2,0,1])
		green = color_own([0,1,0,1])
		red = color_own([1,0,0,1])
		orange = color_own([1,1,0,1]) 

	else: 
		black = color_own([0.,0.,0.,1])
		grey = color_own([.5,.5,0.5,1])
		cyan = color_own([0,1,1,1])
		blue = color_own([0,0,1,1])
		lime = color_own([0.6,1.2,0,1])
		green = color_own([0,1,0,1])
		red = color_own([1,0,0,1])
		orange = color_own([1,1,0,1]) 

	t = np.array([12, 35, 41, 61, 73, 89])
	scandir = np.array([0.1, 0.7, 0.4, 0.8 , 0.2, 0.1])*np.pi	

	#Position_lens
	x1 = np.linspace(1,13,100) + 0.3 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100)) 
	y1 = np.linspace(1,3,100) + 0.3* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	

	#unlensed Position_source
	x2 = np.linspace(5,9,100)# +  0.03 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y2 = np.linspace(7,4.5,100)# + 0.03* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	


	d = np.sqrt((x1-x2)**2 + (y1-y2)**2)
	TE = 2
	X2 = x2 + (x2-x1) * TE/(d**2 +2*TE)
	Y2 = y2 + (y2-y1) * TE/(d**2 +2*TE)	

	dX2 = x1-X2
	dY2 = y1-Y2
	dx2 = x1-x2
	dy2 = y1-y2	

	fig = plt.figure(figsize= (12,8))
	ax = plt.subplot(111)
	ax.axis('equal')
	ax.axis('off')	
	#---------------------------------------------------------------
	#axis 
	plt.xlim([-0.5,14])	
	xr = 12
	yr = 7
	plt.text(xr-0.8,0,'RA $\cdot$ cos(Dec)',verticalalignment = 'center',fontsize = 25)
	plt.text(0,yr + 0.25,'Dec',fontsize = 25, horizontalalignment = 'center', rotation = 90)	

	plt.arrow(-0.025,0,xr-1,0,width = 0.05,overhang = 0.5,head_width = 0.5, head_length = 0.5,color= black, zorder = 100,length_includes_head=True)
	plt.arrow(0,-0.025,0,yr-0.5,width = 0.05,overhang = 0.5,head_width = 0.5,head_length = 0.5,color= black, zorder = 100,length_includes_head=True)	
	
	plt.text(2,1.5,'Lens',color = grey, fontsize = 25, horizontalalignment = 'center', rotation = 0, weight = 'bold')	
	#---------------------------------------------------------------
	# Motion source 
	plt.plot(x1,y1, color = grey, linewidth = 7)
	fig.savefig(imagepath + string + '_1.png', format = 'png')	
	plt.text(4,7.5,'Source',color = blue,fontsize = 25, horizontalalignment = 'center', rotation = 0,weight = 'bold')	

	plt.plot(x2,y2, color = cyan, linestyle= '--',dashes=(10, 5), linewidth = 3, zorder = -1)
	fig.savefig(imagepath + string + '_2.png', format = 'png')	
	plt.plot(X2,Y2, color = blue, linewidth = 3)
	for i in range(len(t)): 
		plt.plot([x2[t[i]],X2[t[i]]],[y2[t[i]],Y2[t[i]]],':',color = black)	
	fig.savefig(imagepath + string + '_3.png', format = 'png')



	delta = 0.05
	for i in range(len(t)): 
		xm1 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + x1[t[i]] 
		ym1 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + y1[t[i]]
		xm2 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + X2[t[i]]
		ym2 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + Y2[t[i]]
		plt.plot(xm2,ym2, color = black, linewidth = 1,zorder = 1)
		plt.plot(xm1,ym1, color = black, linewidth = 1,zorder = 1)
	fig.savefig(imagepath + string + '_4.png', format = 'png')

	for i in range(len(t)): 
		xm1 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + x1[t[i]] 
		ym1 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + y1[t[i]]
		xm2 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + X2[t[i]]
		ym2 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + Y2[t[i]]

		dsc = ((dx2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dy2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]
		dSC = ((dX2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dY2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]	

		ttX2 = np.array([0,-dSC[1]/2,dSC[1]/2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],X2[t[i]],X2[t[i]]])
		ttY2 = np.array([0,-dSC[1]/2,dSC[1]/2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],Y2[t[i]],Y2[t[i]]])
		ttx2 = np.array([0,-dsc[1]/2-0.2,dsc[1]/2-0.2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],x2[t[i]],x2[t[i]]])
		tty2 = np.array([0,-dsc[1]/2-0.2,dsc[1]/2-0.2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],y2[t[i]],y2[t[i]]])	

		if i % 2  == 0:
			plt.arrow(ttx2[2],tty2[2], 0.0001*(ttx2[2]-ttx2[1]),0.0001*(tty2[2]-tty2[1]),color= red, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttx2[1],tty2[1], 0.0001*(ttx2[1]-ttx2[2]),0.0001*(tty2[1]-tty2[2]),color= red, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[2],ttY2[2], 0.0001*(ttX2[2]-ttX2[1]),0.0001*(ttY2[2]-ttY2[1]),color= orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[1],ttY2[1], 0.0001*(ttX2[1]-ttX2[2]),0.0001*(ttY2[1]-ttY2[2]),color= orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True, zorder = 10)
			plt.plot(ttx2[0:2],tty2[0:2],color = black, linestyle= ':')
			plt.plot(ttX2[0:2],ttY2[0:2],color = black, linestyle= ':')
			plt.plot(ttx2[1:3],tty2[1:3],color = red,linewidth = 3 ,  linestyle= '--',dashes=(10, 10))
			plt.plot(ttX2[1:3],ttY2[1:3],color = orange, linewidth = 3,linestyle= '-')
			plt.plot(ttx2[2:],tty2[2:],color = black, linestyle= ':')
			plt.plot(ttX2[2:],ttY2[2:],color = black, linestyle= ':')
		#if i ==0 :



	fig.savefig(imagepath + string + '_5.png', format = 'png')	
	print('Create Image: '+ imagepath+ string + '.png')
	
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	
	plt.close(fig)

def Image_Illustration_Multi(string = 'Illustration_Multi'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
		grey = color_own([.5,.5,0.5,1])
		cyan = color_own([0,1,1,1])
		blue = color_own([0,0,1,1])
		lime = color_own([0.6,1.2,0,1])
		green = color_own([0,1,0,1])
		red = color_own([1,0,0,1])
		orange = color_own([1,.5,0,1]) 
	else:
		black = color_own([0.,0.,0.,1])
		grey = color_own([.5,.5,0.5,1])
		cyan = color_own([0,1,1,1])
		blue = color_own([0,0,1,1])
		lime = color_own([0.6,1.2,0,1])
		green = color_own([0,1,0,1])
		red = color_own([1,0,0,1])
		orange = color_own([1,1,0,1]) 
	t = np.array([12, 35, 41, 61, 73, 89])
	scandir = np.array([0.1, 0.7, 0.4, 0.8 , 0.2, 0.1])*np.pi	

	x1 = np.linspace(1,13,100) + 0.3 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y1 = np.linspace(1,3,100) + 0.3* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	

	x2 = np.linspace(3,7,100)# +  0.03 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y2 = np.linspace(7,4.5,100)# + 0.03* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))	

	x3 = np.linspace(12,10,100)# +  0.03 * np.sin(np.linspace(np.pi/4,12*np.pi/4,100))
	y3 = np.linspace(8,6,100)# + 0.03* np.cos(np.linspace(np.pi/4,12*np.pi/4,100))

	d2 = np.sqrt((x1-x2)**2 + (y1-y2)**2)
	d3 = np.sqrt((x1-x3)**2 + (y1-y3)**2)

	TE = 1.5
	X2 = x2 + (x2-x1) * TE/(d2**2 +2*TE)
	Y2 = y2 + (y2-y1) * TE/(d2**2 +2*TE)	
	X3 = x3 + (x3-x1) * TE/(d3**2 +2*TE)
	Y3 = y3 + (y3-y1) * TE/(d3**2 +2*TE)	

	dX2 = x1-X2
	dY2 = y1-Y2
	dx2 = x1-x2
	dy2 = y1-y2	

	dX3 = x1-X3
	dY3 = y1-Y3
	dx3 = x1-x3
	dy3 = y1-y3	

	fig = plt.figure(figsize= (12,8.8))
	plt.subplots_adjust(
		top=0.95,
		bottom=0.05,
		left=0.05,
		right=0.95,
		hspace=0.0,
		wspace=0.2)	
	ax = plt.subplot(111)
	ax.axis('equal')
	ax.axis('off')	
	
	for i in range(len(t)): 
		delta = 0.05
		xm1 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + x1[t[i]] 
		ym1 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + y1[t[i]]
		xm2 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + X2[t[i]]
		ym2 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + Y2[t[i]]	
		xm3 =np.array([-1,1,1,-1,-1]) * np.cos(scandir[i]) + delta * np.array([-1,-1,1,1,-1]) * np.sin(scandir[i]) + X3[t[i]]
		ym3 =np.array([-1,1,1,-1,-1]) * np.sin(scandir[i]) - delta * np.array([-1,-1,1,1,-1]) * np.cos(scandir[i]) + Y3[t[i]]	

		dSC3 = ((dX3[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dY3[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]
		dSC2 = ((dX2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dY2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]	

		dsc3 = ((dx3[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dy3[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]
		dsc2 = ((dx2[t[i]]).reshape(-1,1)*[np.sin(scandir[i]),np.cos(scandir[i])] \
			+ (dy2[t[i]]).reshape(-1,1) *[-np.cos(scandir[i]),np.sin(scandir[i])])[0]	

		ttX2 = np.array([0,-dSC2[1]/2,dSC2[1]/2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],X2[t[i]],X2[t[i]]])
		ttY2 = np.array([0,-dSC2[1]/2,dSC2[1]/2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],Y2[t[i]],Y2[t[i]]])
		ttX3 = np.array([0,-dSC3[1]/2,dSC3[1]/2,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],X3[t[i]],X3[t[i]]])
		ttY3 = np.array([0,-dSC3[1]/2,dSC3[1]/2,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],Y3[t[i]],Y3[t[i]]])	
		ttx2 = np.array([0,-dsc2[1]/2-0.3,dsc2[1]/2-0.3,0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],x2[t[i]],x2[t[i]]])
		tty2 = np.array([0,-dsc2[1]/2-0.3,dsc2[1]/2-0.3,0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],y2[t[i]],y2[t[i]]])
		
		if i == 3: off = [-0.4,-0.2]
		else: off = [0,-0.2]
		ttx3 = np.array([0,-dsc3[1]/2+off[1],dsc3[1]/2+off[1],0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],x3[t[i]],x3[t[i]]])
		tty3 = np.array([0,-dsc3[1]/2+off[1],dsc3[1]/2+off[1],0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],y3[t[i]],y3[t[i]]])	
		ttX3 = np.array([0,-dSC3[1]/2+off[0],dSC3[1]/2+off[0],0]) * np.cos(scandir[i]) + ([x1[t[i]],x1[t[i]],X3[t[i]],X3[t[i]]])
		ttY3 = np.array([0,-dSC3[1]/2+off[0],dSC3[1]/2+off[0],0])  * np.sin(scandir[i]) +([y1[t[i]],y1[t[i]],Y3[t[i]],Y3[t[i]]])	

		'''
		if i % 2  == 0:
			plt.arrow(ttX2[2],ttY2[2], 0.0001*(ttX2[2]-ttX2[1]),0.0001*(ttY2[2]-ttY2[1]),color =  color_own([0,0,1,1]),head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[1],ttY2[1], 0.0001*(ttX2[1]-ttX2[2]),0.0001*(ttY2[1]-ttY2[2]),color =  color_own([0,0,1,1]),head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)

			plt.plot(ttX2[0:2],ttY2[0:2],color = black, linestyle= ':')
			plt.plot(ttX2[1:3],ttY2[1:3],color = [158/200,1/200,66/200, 1], linewidth = 3,linestyle= '-')
			plt.plot(ttX2[2:],ttY2[2:],color = black, linestyle= ':')
		'''
		if i% 2 == 0: 
			plt.plot(xm2,ym2, color = black, linewidth = 1,zorder = 1)
			plt.plot(xm1,ym1, color = black, linewidth = 1,zorder = 1)

			plt.arrow(ttX2[2],ttY2[2], 0.0001*(ttX2[2]-ttX2[1]),0.0001*(ttY2[2]-ttY2[1]),color=red ,head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX2[1],ttY2[1], 0.0001*(ttX2[1]-ttX2[2]),0.0001*(ttY2[1]-ttY2[2]),color=red ,head_width=0.2,\
				overhang = 0.5, length_includes_head=True, zorder = 10)


			plt.plot(ttX2[0:2],ttY2[0:2],color = black, linestyle= ':')
			plt.plot(ttX2[1:3],ttY2[1:3],color = red, linewidth = 3)
			plt.plot(ttX2[2:],ttY2[2:],color = black, linestyle= ':')
			
			plt.plot(ttx2[0:2],tty2[0:2],color = black, linestyle= ':')
			plt.plot(ttx2[1:3],tty2[1:3],color =  orange, linewidth = 3,linestyle= '-', dashes=(10, 2))
			plt.plot(ttx2[2:],tty2[2:],color = black, linestyle= ':')
			plt.arrow(ttx2[2],tty2[2], 0.0001*(ttx2[2]-ttx2[1]),0.0001*(tty2[2]-tty2[1]),color =  orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttx2[1],tty2[1], 0.0001*(ttx2[1]-ttx2[2]),0.0001*(tty2[1]-tty2[2]),color = orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True, zorder = 10)
			if i >=3:
				plt.plot(ttX3[0:2],ttY3[0:2],color = black, linestyle= ':')
				plt.plot(ttX3[1:3],ttY3[1:3],color = red, linewidth = 3)
				plt.plot(ttX3[2:],ttY3[2:],color = black, linestyle= ':')
				plt.arrow(ttX3[2],ttY3[2], 0.0001*(ttX3[2]-ttX3[1]),0.0001*(ttY3[2]-ttY3[1]),color=red, head_width=0.2,\
					overhang = 0.5, length_includes_head=True ,zorder = 10)
				plt.arrow(ttX3[1],ttY3[1], 0.0001*(ttX3[1]-ttX3[2]),0.0001*(ttY3[1]-ttY3[2]),color= red, head_width=0.2,\
					overhang = 0.5, length_includes_head=True ,zorder = 10)
				plt.plot(xm3,ym3, color = black, linewidth = 1,zorder = 1)

				plt.plot(ttx3[0:2],tty3[0:2],color = black, linestyle= ':')
				plt.plot(ttx3[1:3],tty3[1:3],color =  orange, linewidth = 3,linestyle= '-', dashes=(10, 2))
				plt.plot(ttx3[2:],tty3[2:],color = black, linestyle= ':')
				plt.arrow(ttx3[2],tty3[2], 0.0001*(ttx3[2]-ttx3[1]),0.0001*(tty3[2]-tty3[1]),color =  orange, head_width=0.2,\
					overhang = 0.5, length_includes_head=True ,zorder = 10)
				plt.arrow(ttx3[1],tty3[1], 0.0001*(ttx3[1]-ttx3[2]),0.0001*(tty3[1]-tty3[2]),color =  orange, head_width=0.2,\
					overhang = 0.5, length_includes_head=True, zorder = 10)
			'''
			else:
				plt.plot(xm3,ym3, color = 'grey', linewidth = 2, zorder = -1)
			'''

		elif i >=3: 
			plt.plot(xm3,ym3, color = black, linewidth = 1,zorder = 1)
			plt.plot(xm1,ym1, color = black, linewidth = 1,zorder = 1) 

			#plt.plot(xm2,ym2, color = 'grey', linewidth = 2, zorder = -1)

			plt.plot(ttX3[0:2],ttY3[0:2],color = black, linestyle= ':')
			plt.plot(ttX3[1:3],ttY3[1:3],color = red, linewidth = 3 )
			plt.plot(ttX3[2:],ttY3[2:],color = black, linestyle= ':')
		
			plt.arrow(ttX3[2],ttY3[2], 0.0001*(ttX3[2]-ttX3[1]),0.0001*(ttY3[2]-ttY3[1]),color= red, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttX3[1],ttY3[1], 0.0001*(ttX3[1]-ttX3[2]),0.0001*(ttY3[1]-ttY3[2]),color= red, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)

			plt.plot(ttx3[0:2],tty3[0:2],color = black, linestyle= ':')
			plt.plot(ttx3[1:3],tty3[1:3],color =  orange, linewidth = 3,linestyle= '-', dashes=(10, 2))
			plt.plot(ttx3[2:],tty3[2:],color = black, linestyle= ':')
			plt.arrow(ttx3[2],tty3[2], 0.0001*(ttx3[2]-ttx3[1]),0.0001*(tty3[2]-tty3[1]),color= orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True ,zorder = 10)
			plt.arrow(ttx3[1],tty3[1], 0.0001*(ttx3[1]-ttx3[2]),0.0001*(tty3[1]-tty3[2]),color= orange, head_width=0.2,\
				overhang = 0.5, length_includes_head=True, zorder = 10)


		'''
		else:  
			plt.plot(xm3,ym3, color = 'grey', linewidth = 2, zorder = -1)
			plt.plot(xm2,ym2, color = 'grey', linewidth = 2, zorder = -1)
			plt.plot(xm1,ym1, color = 'grey', linewidth = 2, zorder = -1)
		'''
		#if i ==0 :
	plt.plot(x1,y1, color = grey, linewidth = 7)
	plt.plot(x2,y2, color = cyan,  linestyle= '--',dashes=(6, 2), linewidth = 3, zorder = -1)
	plt.plot(X2,Y2, color = blue, linewidth = 3)	
	plt.plot(x3,y3, color = lime,  linestyle= '--',dashes=(6, 2), linewidth = 3, zorder = -1)
	plt.plot(X3,Y3, color = green, linewidth = 3)	

	plt.xlim([-0.2,13.5])	

	xr = 12
	yr = 7
	plt.text(xr-0.8,0,'RA $\cdot$ cos(Dec)',verticalalignment = 'center',fontsize = 25)
	plt.text(0,yr + 0.25,'Dec',fontsize = 25, horizontalalignment = 'center', rotation = 90)	
	plt.text(2,1.5,'Lens',color = grey,fontsize = 25, horizontalalignment = 'center', rotation = 0, weight = 'bold')	
	plt.text(4,7.5,'Star 1',color = blue, fontsize = 25, horizontalalignment = 'center', rotation = 0,weight = 'bold')	
	plt.text(11,8,'Star 2',color = green, fontsize = 25, horizontalalignment = 'center', rotation = 0,weight = 'bold')	

	plt.arrow(-0.025,0,xr-1,0,width = 0.05,overhang = 0.5,head_width = 0.5, head_length = 0.5,color= black, zorder = 100,length_includes_head=True)
	plt.arrow(0,-0.025,0,yr-0.5,width = 0.05,overhang = 0.5,head_width = 0.5,head_length = 0.5,color= black, zorder = 100,length_includes_head=True)	

	fig.savefig(imagepath + string + '.png', format = 'png')	
	print('Create Image: '+ imagepath+ string + '.png')
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	plt.close(fig)

def Image_compare_micro(string = 'aml_vs_pml.png'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark:
		string = 'dark_'+string 
		c1 =color_own([1,1,0,1])
		c2 =color_own([1,0,1,1])
		c3 =color_own([0.6,1.2,0,1])
		c4 =color_own([0,1,1,1])
	else:
		c1 =color_own([0,0,1,1])
		c2 =color_own([1,0,0,1])
		c3 =color_own([0.5,1,0,1])
		c4 =color_own([1,0,1,1])

	rc('xtick', labelsize=24) 
	rc('ytick', labelsize=24) 
	Separation_E1 = 20
	Separation_E2 = 2
	Separation_E3 = 20
	xx1 = np.array(range(-1000,1000,1))
	xx2 = xx1
	xx3 = xx1
	yy1 = xx1*0+10
	yy2 = xx1*0+1
	yy3 = xx1*0+200
	
	uu1 = np.sqrt(xx1*xx1+yy1*yy1)/Separation_E1
	uu2 = np.sqrt(xx2*xx2+yy2*yy2)/Separation_E2
	uu3 = np.sqrt(xx3*xx3+yy3*yy3)/Separation_E3
	
	dSeparation1 = uu1/(uu1*uu1 + 2)*Separation_E1
	dSeparation2 = uu2/(uu2*uu2 + 2)*Separation_E2
	dSeparation3 = uu3/(uu3*uu3 + 2)*Separation_E3
	
	A1 = (uu1*uu1+2)/(uu1*np.sqrt(uu1*uu1+4))
	A2 = (uu2*uu2+2)/(uu2*np.sqrt(uu2*uu2+4))
	A3 = (uu3*uu3+2)/(uu3*np.sqrt(uu3*uu3+4))
	
	dm1 = 2.5*np.log10(A1)
	dm2 = 2.5*np.log10(A2)
	dm3 = 2.5*np.log10(A3)
	
	xx1 = xx1/250
	xx2 = xx2/250
	xx3 = xx3/250
	figure = plt.figure(figsize = (12,12))
	plt.subplots_adjust(hspace=0.1)	

	ax1 = plt.subplot2grid((2,1), (0, 0), rowspan=1)
	plt.plot(xx1,dSeparation1, color = c1,linewidth = 4)
	line1, = plt.plot(xx2,dSeparation2,color = c3 ,linewidth = 4)
	plt.plot(xx3,dSeparation3,linestyle = '--',color = c4,linewidth = 4)
	line1.set_dashes([10, 2, 10, 2])	

	plt.yticks([1,2,3,4,5,6,7,8],['1.0 ','2.0 ','3.0 ','4.0 ','5.0 ','6.0 ','7.0 ','8.0 '])
	plt.ylim([0,8])
	plt.ylabel('Shift [mas]',fontsize = 30)
	plt.plot(xx1,xx1*0+0.5, color=c2,linewidth = 3)
	ax1.tick_params(length=6, width=2)
	ax1.tick_params(which='minor', length=4,width=2)
	xticklabels1 = ax1.get_xticklabels()
	plt.setp(xticklabels1, visible=False)
	ax2 = plt.subplot2grid((2,1), (1, 0), rowspan=1,sharex=ax1)	
	plt.semilogy(xx1,dm1,color = c1,linewidth = 4)
	line1, = plt.semilogy(xx2,dm2,color = c3,linewidth = 4)
	plt.semilogy(xx3,dm3,linestyle = '--',color = c4,linewidth = 4)
	line1.set_dashes([10, 2, 10, 2])
	plt.semilogy(xx1,xx1*0+0.001, color= c2  ,linewidth = 3)
	plt.ylim([0.0001,1])
	plt.ylabel('Magnification  [mag]',fontsize = 30)
	ax2.tick_params(length=6, width=2)
	ax2.tick_params(which='minor', length=4,width=2)	
	plt.xlabel('$\Delta$ T [yr]',fontsize = 30)
	figure.savefig(imagepath+ string,format = 'png')
	print('Create Image: '+ imagepath+ string + '.png')

	plt.close(figure)


def Image_astroshift(string = 'astroshift'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	fig,ax = plt.subplots(figsize = [12,12]	)
	x = np.array([-100+0.1*t for t in range(2000)])
	ThetaE = 12.75
	umin = 0.75
	ThetaE = 12.75
	fls  = 10
	siz = 10

	y =  np.array([ThetaE*umin for i in x])
	ax.set_xlim(-30,30)
	ax.set_ylim(-30,30)	

	ux =np.array([x[i]/ThetaE for i in range(len(x))])
	uy =np.array([y[i]/ThetaE for i in range(len(x))])
	u = np.array([np.sqrt(np.power(ux[i],2)+np.power(uy[i],2)) for i in range(len(x))])	
	theta_px = - (np.sqrt(u*u+4) - u)/(2 * u) * x
	theta_py = - (np.sqrt(u*u+4) - u)/(2 * u) * y
	theta_mx = - (-np.sqrt(u*u+4) - u)/(2 * u) * x
	theta_my = - (-np.sqrt(u*u+4) - u)/(2 * u) * y


	

	delta_thetax =  -x/(np.power(u,2)+2)
	delta_thetay =  -y/(np.power(u,2)+2)

	plt.plot(x,y,'r',linewidth = 2)
	plt.plot(theta_px,theta_py,'b--',linewidth = 2)
	plt.plot(theta_mx,theta_my,'b--',linewidth = 2)
	plt.plot(delta_thetax,delta_thetay,color='green',linewidth = 2)
	xx = [ThetaE * np.sin(2*np.pi*i/1000)+6 for i in range(1000)] 
	yy = [ThetaE * np.cos(2*np.pi*i/1000)+ThetaE*umin  for i in range(1000)]	
	step = 30
	x = x[10::step]
	y = y[10::step]
	theta_px = theta_px[10::step]
	theta_py = theta_py[10::step]
	theta_mx = theta_mx[10::step]
	theta_my = theta_my[10::step]
	delta_thetax = delta_thetax[10::step]
	delta_thetay = delta_thetay[10::step]

	plt.plot(x,y,'rx', markersize = siz, markeredgewidth = 3)
	plt.plot(theta_px,theta_py,'bs', markersize = siz)
	plt.plot(theta_mx,theta_my,'bs', markersize = siz)
	plt.plot(delta_thetax,delta_thetay,'o', color='green')
	n = np.where(np.array(x) == 6.0)[0][0]
	plt.plot(x[n],y[n],'rx',markersize = siz * 1.5,markeredgewidth = 4)

	plt.plot(theta_px[n],theta_py[n],'bs',markersize = siz * 1.5)
	plt.plot(theta_mx[n],theta_my[n],'bs',markersize = siz * 1.5)
	plt.plot(delta_thetax[n],delta_thetay[n],'o', color='green',markersize = siz * 1.5)

	plt.plot(xx,yy,color='black')
	plt.plot([-60,60],[-10*ThetaE*umin ,10 * ThetaE*umin ], '--', color = 'black', dashes=[10,5])
	plt.xlabel(r'$\delta\theta_{x} [\mathrm{mas}]$', fontsize = 30)
	plt.ylabel(r'$\delta\theta_{y} [\mathrm{mas}]$', fontsize = 30)

	plt.arrow(13.5,21,2,-1,width = 0.1,head_width = 0.70,head_length = 1,color= 'black', zorder = 100)
	plt.arrow(6.75,10.75,2,0,width = 0.1,head_width = 0.70,head_length = 1,color= 'red', zorder = 100)
	plt.arrow(11,18,2,-0.6,width = 0.1,head_width = 0.70,head_length = 1,color= 'blue', zorder = 100)
	plt.arrow(-1,-5,-2,0.7,width = 0.1,head_width = 0.7,head_length = 1,color= 'green', zorder = 100)
	plt.arrow(-4,-9,-2,1,width = 0.1,head_width = 0.7,head_length = 1,color= 'blue', zorder = 100)

	#plt.scatter([0],[0],s = 100, c='k', zorder = 1000)
	plt.plot([0],[0],'k*',markersize = siz * 2)

	plt.text(-10,20, '$-$',color = 'blue', fontsize = 50)
	plt.text(7,-10, '$+$',color = 'blue',fontsize = 50)	
	plt.xticks( fontsize = 25)
	plt.yticks( fontsize = 25)
	fig.savefig(imagepath + string + '.png')
	if paperpath is not None: fig.savefig(paperpath+ string + '.png')	
	print('Create Image: '+ imagepath+ string + '.png')

def create_all_Image():
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	# code to create all images in default mode at once
	for j in range(1,-1,-1):
		dark(j)
		for i in __all__: #loop over all functions in this code 
			if 'Image_' in i:
				exec(i+'()')

				
if __name__ == '__main__': create_all_Image()

