import matplotlib.pyplot as plt
import numpy as np
__all__ = ['color_own']


def color_own(length, rng = [0,1], colormap = plt.cm.gist_rainbow, ):
	m = 1
	w = np.array([245,245,245,0]).reshape(-1,4)
	m = np.array([145,22,104,0]).reshape(-1,4)
	c = np.array([0,176,230,0]).reshape(-1,4)
	y = np.array([250,194,0,0]).reshape(-1,4)
	r = np.array([225,0,35,0]).reshape(-1,4)
	g = np.array([1,143,53,0]).reshape(-1,4)
	b = np.array([0,66,148,0]).reshape(-1,4)
	alpha = [0,0,0,255]
	if isinstance(length, int):
		color_give = False
		c1 = colormap(np.linspace(*rng,length)).reshape(-1,4)
	else: 
		color_give = True	
		c1 = np.array(length).reshape(-1,4)
	cmin = np.argmin(c1[:,:3],axis = 1)
	cmax = np.argmax(c1[:,:3],axis = 1)
	case1 = np.where((cmin == 0))
	case2 = np.where((cmin == 1))
	case3 = np.where((cmin == 2))
	case4 = np.where((cmax == 0))
	case5 = np.where((cmax == 1))
	case6 = np.where((cmax == 2))
	wms = np.sort(c1[:,:3],axis = 1)
	white  =  wms[:,0].reshape(-1,1)
	wms = wms -white
	mix = wms[:,1].reshape(-1,1)
	wms = wms -mix
	col = wms[:,2].reshape(-1,1)
	color = w * white
	color[case1] += c * mix[case1]
	color[case2] += m * mix[case2]
	color[case3] +=	y * mix[case3]
	color[case4] += r * col[case4]
	color[case5] += g * col[case5]
	color[case6] += b * col[case6]
	color += alpha
	if color_give: 
		if len(color) == len(length):
			return color/255 
		else: return color.reshape(4)/255
	return color/255