import numpy as np
import glob,os
import scipy.io
#~ from scipy.io.matlab import mio
from scipy import io
from scipy import ndimage as nd, misc
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy.ma as ma
from scipy import interpolate
import threading

def get_steps(data_dir):
	#~ This function gives back an array with the steps.
	list_steps = np.sort(np.array(glob.glob('%s/*mat' % data_dir)),axis=0)
	steps = np.zeros([np.shape(list_steps)[0]],dtype=int)
	for i,step in enumerate(list_steps):
		#~ steps[i]=int(step[-11:-6])
		steps[i]=int(step.split('/')[-1].split('-')[-1][:-6])
	step_0,step_f = get_start_end_steps(data_dir)
	
	steps = steps[(steps<=step_f)*(steps>=step_0)]
	
	return steps

def get_basename(data_dir):
	#~ This function gives back a string with the root name of every correlation step.
	#~ Ex: 'tensile_test_steel01-rup_00001.mat'  ->  'tensile_test_steel01-rup'
	#~ print np.sort(np.array(glob.glob('%s/*mat' % data_dir)),axis=0)[0].split('/')[-1][:-12]
	
	return(np.sort(np.array(glob.glob('%s/*mat' % data_dir)),axis=0)[0].split('/')[-1].split('-')[0])

def get_sens(data_dir):
	file_path = data_dir + get_basename(data_dir) + '.dat' 
	
	if os.path.isfile(file_path) == False:
		with open(file_path, 'w') as datafile:
			datafile.write('# 0-x	1-y	2-angle	3-pos	4-sens	5-dx	6-dy	7-step_0	8-step_f\n')
			datafile.write('1.	1.	0.	0.	1.	0.	0.	0.	1000.\n')
	return np.loadtxt(file_path)[4]	
	
def get_field_params(data_dir):
	file_path = data_dir + get_basename(data_dir) + '.dat' 
	
	if os.path.isfile(file_path) == False:
		with open(file_path, 'w') as datafile:
			datafile.write('# 0-x	1-y	2-angle	3-pos	4-sens	5-dx	6-dy	7-step_0	8-step_f\n')
			datafile.write('1.	1.	0.	0.	1.	0.	0.	0.	1000.\n')

	data	= 	np.loadtxt(file_path)
	x		=	data[0]
	y		=	data[1]
	angle	=	data[2]
	pos 	=	data[3]
	sens 	=	data[4]
	dx	 	=	data[5]
	dy	 	=	data[6]
	step_0	=	data[7]
	step_f 	=	data[8]
	
	#~ return x,y,angle,pos,sens,dx,dy
	return x,y,angle,pos,sens,dx,dy,step_0,step_f

def get_field(data_dir,field,step):
	#~ The function io.savemat handles the dict names oddly. This does not allow us to store in the disk the corrected data.
	#~ Thus every time data is read, it needs to be corrected.
	#~ x,y,angle,pos,sens,dx,dy = get_field_params(data_dir)
	if data_dir[-1] <> "/":
		data_dir += '/'
	
	x,y,angle,pos,sens,dx,dy,step_0,step_f = get_field_params(data_dir)
	basename = get_basename(data_dir)
	
	#~ print('Opening field: %s' % field)

	if step > get_steps(data_dir)[-1]:
		step = str(int(step_f))
		#~ print step
	
	nz = len(glob.glob(data_dir+'*00_0.mat')[0].split('/')[-1].split('-')[-1][:-6])
	file_name = data_dir + basename + '-' + str(step).zfill(nz) +'_0.mat'
	
	Z0 = np.array((io.loadmat(file_name))[field],dtype=np.float64)
	sigma = np.array((io.loadmat(file_name))['sigma'],dtype=np.float64)
	mask = sigma<=0
	data= ma.masked_array(Z0, mask=mask)
	
	'''
	sig = np.array((io.loadmat(file_name))['sigma'],dtype=np.float64)
	#~ data= Z0
	plt.imshow(sig,interpolation='nearest')
	plt.show()
	
	x = np.array((io.loadmat(file_name))['X'],dtype=np.float64)
	
	y = np.array((io.loadmat(file_name))['Y'],dtype=np.float64)
	
	plt.plot(y[::,10],y[::,10])
	plt.plot(ma.masked_array(y, mask=mask)[::,10],ma.masked_array(y, mask=mask)[::,10],'o')
	plt.show()
	

	#~ plt.plot(y[::,10],sigma[::,10])
	plt.plot(ma.masked_array(y, mask=mask)[::,10],ma.masked_array(sigma, mask=mask)[::,10],'o')
	plt.plot(y[::,10],sigma[::,10],'o')
	plt.show()
	
	Z0 = np.array((io.loadmat(file_name))[field],dtype=np.float64)
	plt.subplot(141)
	#~ plt.imshow(x,interpolation='nearest')
	plt.imshow(ma.masked_array(x, mask=mask),interpolation='nearest')
	plt.subplot(142)
	plt.imshow(y,interpolation='nearest')
	plt.subplot(143)
	plt.imshow(sigma,interpolation='nearest')
	plt.subplot(144)
	plt.imshow(mask,interpolation='nearest')
	plt.show()
	'''
	#~ data = ndimage.rotate(Z0, angle, axes=(1, 0), reshape=False)
	if field=='X' or field=='x_c': 
		data = data-dx
		
	if field=='Y' or field=='y_c': 
		data+=-dy
	return data[x:-x:,y:-y:]

def check_correlations(data_dir):
	result = 'This correlations looks fine to me! =)'
	# This function checks if every step of a given correlation has been entirely correlated. Sometimes correlations fails and there are nan holes in the output field.
	steps=get_steps(data_dir)
	for step in steps:
		#~ X,Y,Z = get_field(data_dir,'X',step)
		eyy = get_field(data_dir,'eyy',step)
		if np.mean(np.isnan(eyy))!=0:
			print 'Step %s has not been completely correlated.  ' % step
			result = 'This correlations needs to be checked! =/'
	print result
	return

def get_mesh(data_dir):
	step = get_steps(data_dir)[0]
	
	file_name = data_dir +get_basename(data_dir)+'-' + str(step).zfill(5) +'_0.mat'

	X = np.array((io.loadmat(file_name))['X'])[0,::]
	Y = np.array((io.loadmat(file_name))['Y'])[::,0]
	X, Y = np.meshgrid(X, Y)
	return X,Y

def get_mesh_size(data_dir):
	steps = get_steps(data_dir)
	step=steps[0]
	X = get_field(data_dir,'X',steps[0])
	Ty,Tx = np.shape(X)
	return Tx,Ty
	
def set_corr_window():
	n = 0
	Active = True
	while Active:
		os.system('clear')
		if n==0:
			print 'Enter coordinates of the first point in the reference axis (separated by coma):'
			p2=raw_input()
			
			print 'Enter coordinates of the second point in the reference axis (separated by coma):'
			p3=raw_input()
			
			x2 = int(p2.split(',')[0])
			y2 = int(p2.split(',')[1])

			x3 = int(p3.split(',')[0])
			y3 = int(p3.split(',')[1])
			
			P2=np.array([x2,y2])
			P3=np.array([x3,y3])
			
			v=P3-P2
			u=np.array([v[1],-v[0]])
			n=1
		
		print 'Enter coordinates of the starting point (separated by coma):'
		p0=raw_input()
		
		x0 = int(p0.split(',')[0])
		y0 = int(p0.split(',')[1])
		
		P0=np.array([x0,y0])
		#~ print 'P0 = %s' % P0
		#~ print 'P2 = %s' % P0
		#~ print 'P3 = %s' % P3
				
		dx=int(raw_input('Width? '))
		i=u/np.linalg.norm(u)*dx
		dy=int(raw_input('Height? '))
		j=v/np.linalg.norm(v)*dy
		
		NW = (P0).astype(int)
		SW = (P0+j).astype(int)
		SE = (P0+i+j).astype(int)
		NE = (P0+i).astype(int)
		
		print 'Coordinates of the window:'
		print 'NW = %s' % NW
		print 'SW = %s' % SW
		print 'SE = %s' % SE
		print 'NE = %s' % NE
		
		Active = True if (raw_input('Do you want to continue? (y/n)')=='y') else False
	if (raw_input('Do you want to save this settings? (y/n)')=='y'):
		from datetime import datetime, date, time
		weld_name = raw_input('Weld name? ')
		datafile = '00-correlations_windows/'+datetime.now().isoformat() + '--' + weld_name + '.dat'
		with open(datafile, 'a') as datafile:
			datafile.write('Weld: %s \n\n' % weld_name)
			datafile.write('P0 = %s \n' % P0)
			datafile.write('P2 = %s \n' % P2)
			datafile.write('P3 = %s \n\n' % P3)
			datafile.write('Window width = %s\n' % dx)
			datafile.write('Window height = %s\n\n' % dy)
			datafile.write('NW = %s\n' % NW)
			datafile.write('SW = %s\n' % SW)
			datafile.write('SE = %s\n' % SE)
			datafile.write('NE = %s\n\n' % NE)
			datafile.write('Axis vectors: \n')
			datafile.write('i = %s \n' % i)
			datafile.write('j = %s \n' % j)
			
		
	return
	
def get_params(analysis_params):	
	i	= int(analysis_params.split('_')[0])
	ss	= int(analysis_params.split('_')[1])
	st	= int(analysis_params.split('_')[2])
	sw	= int(analysis_params.split('_')[3])
	vsg	= int(analysis_params.split('_')[4])
	return i,ss,st,sw,vsg

def get_xy_mm(data_dir):
	steps=get_steps(data_dir)
	x = get_field(data_dir,'X',steps[0])
	y = get_field(data_dir,'Y',steps[0])
	return x[0,::],y[::,0]

def get_extent(data_dir):
	steps=get_steps(data_dir)
	try:
		x = get_field(data_dir,'X',steps[0])
		y = get_field(data_dir,'Y',steps[0])
	except:
		x = get_field(data_dir,'x_c',steps[0])
		y = get_field(data_dir,'y_c',steps[0])
	
	return np.ma.MaskedArray.min(x),np.ma.MaskedArray.max(x),np.ma.MaskedArray.min(y),np.ma.MaskedArray.max(y)
	
def mean_lateral_deformation(data_dir,step):
	y = get_xy_mm(data_dir)[1]
	n = np.shape(y)[0]
	eyy = get_field(data_dir,'eyy',step)
	eyy_mean = np.zeros([n])
	for i in range(n):
		eyy_mean[i] = np.mean(eyy[i,::])
	return y,eyy_mean

def get_value(data_dir,field,step,pt,show=False,verbose=False):
	def pt_inside():	
		contour = np.array([	[xx[np.shape(X)[1]*(np.shape(X)[0]-1)]		,yy[np.shape(X)[1]*(np.shape(X)[0]-1)]	],
								[xx[np.shape(X)[1]*np.shape(X)[0]-1]		,yy[np.shape(X)[1]*np.shape(X)[0]-1]	],
								[xx[np.shape(X)[1]-1]						,yy[np.shape(X)[1]-1]					],
								[xx[0]										,yy[0]									],
								[xx[np.shape(X)[1]*(np.shape(X)[0]-1)]		,yy[np.shape(X)[1]*(np.shape(X)[0]-1)]	]
							])
		p0 = contour[0,::]; pf = contour[1,::]
		L1 = (pf[1]-p0[1])*(pt[0]-p0[0])-(pt[1]-p0[1])*(pf[0]-p0[0])<0
		p0 = contour[1,::]; pf = contour[2,::]
		L2 = (pf[1]-p0[1])*(pt[0]-p0[0])-(pt[1]-p0[1])*(pf[0]-p0[0])<0
		p0 = contour[2,::]; pf = contour[3,::]
		L3 = (pf[1]-p0[1])*(pt[0]-p0[0])-(pt[1]-p0[1])*(pf[0]-p0[0])<0
		p0 = contour[3,::]; pf = contour[0,::]
		L4 = (pf[1]-p0[1])*(pt[0]-p0[0])-(pt[1]-p0[1])*(pf[0]-p0[0])<0	
		return L1*L2*L3*L4

	def get_nodes_around(xx,yy,ff,pt,show):	
		def get_nodes():
			d = np.argmin(np.sqrt((xx-pt[0])**2+(yy-pt[1])**2))
			
			if pt[0]>xx[d] and pt[1]>yy[d]: # first quadrant		
				i1 = d
				i2 = d+1
				i3 = d-np.shape(X)[1]+1
				i4 = d-np.shape(X)[1]
				
			if pt[0]<xx[d] and pt[1]>yy[d]: # second quadrant
				i1 = d-1
				i2 = d
				i3 = d-np.shape(X)[1]
				i4 = d-np.shape(X)[1]-1
				
			if pt[0]<xx[d] and pt[1]<yy[d]: # third quadrant
				i1 = d+np.shape(X)[1]-1
				i2 = d+np.shape(X)[1]
				i3 = d
				i4 = d-1
				
			if pt[0]>xx[d] and pt[1]<yy[d]: # fourth quadrant
				i1 = d+np.shape(X)[1]
				i2 = d+np.shape(X)[1]+1
				i3 = d+1
				i4 = d
			return i1,i2,i3,i4
			
		ele_xy = np.zeros([4,2])
		f = np.zeros([4,1])
		
		
		i1,i2,i3,i4 = get_nodes()
		
		ele_xy = np.array([	[xx[i1],yy[i1]],
							[xx[i2],yy[i2]],
							[xx[i3],yy[i3]],
							[xx[i4],yy[i4]]])
		
		f = np.array([	[ff[i1]],
						[ff[i2]],
						[ff[i3]],
						[ff[i4]] ])
		
		#~ print ele_xy
		#~ print f
		if verbose == True:
			for i in range(4):
				print i+1,ele_xy[i,::],f[i,::]
		return ele_xy,f

	def show_mapping():
		fig = plt.figure()
		plt.title(field)
		plt.plot(xx,yy,'bx')
		plt.imshow(F,extent=[np.min(xx),np.max(xx),np.min(yy),np.max(yy)])
		plt.plot(pt[0],pt[1],'ro')
		try:
			plt.plot(ele_xy[::,0],ele_xy[::,1],'sk')
			plt.text(ele_xy[0,0]-0.05, ele_xy[0,1]-0.05, r'$1$')
			plt.text(ele_xy[1,0]+0.05, ele_xy[1,1]-0.05, r'$2$')
			plt.text(ele_xy[2,0]+0.05, ele_xy[2,1]+0.05, r'$3$')
			plt.text(ele_xy[3,0]-0.05, ele_xy[3,1]+0.05, r'$4$')
		except:
			pass
		plt.axis('equal')
		plt.xlim([pt[0]-0.75,pt[0]+0.75])
		plt.ylim([pt[1]-0.75,pt[1]+0.75])
		plt.colorbar()
		
		plt.show()
		return

	def interpolate_value(ele_xy,ele_z,pt):
		def shape(xi):
			return np.array([[0.25*(1.0-xi[0])*(1.0-xi[1]),0.25*(1.0+xi[0])*(1.0-xi[1]),0.25*(1.0+xi[0])*(1.0+xi[1]),0.25*(1.0-xi[0])*(1.0+xi[1])]])
		def newton(ele,pt):
			xi = np.array([0.0, 0.0])
			for i in range(1,100):
				K = np.dot(dshape(xi),ele)
				b = (pt-np.dot(shape(xi),ele))
				Ki = np.linalg.inv(K)
				dxi = np.dot(b,Ki)
				#~ dxi = K\b;
				if (np.linalg.norm(dxi) < 1.e-8):
					break
			return (xi+dxi)[0]
		def dshape(xi):
			dphi = np.array([	[-0.25*(1.0-xi[1]),  0.25*(1.0-xi[1]),  0.25*(1.0+xi[1]), -0.25*(1.0+xi[1])],
								[-0.25*(1.0-xi[0]), -0.25*(1.0+xi[0]),  0.25*(1.0+xi[0]),  0.25*(1.0-xi[0])]])
			return dphi
		xi_gp = newton(ele_xy,pt);
		return np.dot(shape(xi_gp),ele_z)[0][0]	
	
	
	X = get_field(data_dir,'X',step)
	Y = get_field(data_dir,'Y',step)
	F = get_field(data_dir,field,step)

	xx = X.reshape([np.shape(X)[0]*np.shape(X)[1]])
	yy = Y.reshape([np.shape(Y)[0]*np.shape(Y)[1]])
	ff = F.reshape([np.shape(F)[0]*np.shape(F)[1]])
	
	
	
	if pt_inside():
		ele_xy,f = get_nodes_around(xx,yy,ff,pt,show)
		val = interpolate_value(ele_xy,f,pt)
		if verbose==True:
			print 'Interpolated value : %f.4.5' % val
		if show:
			show_mapping() 
	else:
		print 'This point can not be calculated'
		show_mapping()
		val = np.nan
	
	plt.close()
	return val

def get_close_node(data_dir,pt):
	#~ This function will give the coordinates in pixels of the closest point to the asked one.
	
	step = get_steps(data_dir)[0]
	try:
		X = get_field(data_dir,'X',step)
		Y = get_field(data_dir,'Y',step)
	except:
		X = get_field(data_dir,'x_c',step)
		Y = get_field(data_dir,'y_c',step)

	d = (X-pt[0])**2+(Y-pt[1])**2

	x0,y0  = np.unravel_index(np.nanargmin(d),np.shape(d))
	
	return np.array([x0,y0])

def get_value_beta(data_dir,field,step,pt):
	x0,y0 = get_close_node(data_dir,pt)
	value = get_field(data_dir,field,step)[x0,y0]
	return value
			
def get_value_coord_px(data_dir,field,step,pt):
	x0,y0 = pt
	value = get_field(data_dir,field,step)[x0,y0]
	return value
	
def plot_profile(data_dir,n,label,show=False,xlabel=False,ylim=0.01,color='b',alpha=1):
	steps	= get_steps(data_dir)[:n]
		
	Tx,Ty = get_mesh_size(data_dir)
	
	for i,step in enumerate(steps):
		y 	= get_field(data_dir,'Y',step)[::,Tx//5*4]
		eyy = get_field(data_dir,'eyy',step)[::,Tx//5*4]
		#~ y,eyy = mean_lateral_deformation(data_dir,step)
		
		if i ==0:
			plt.plot(y*get_sens(data_dir),eyy,'%s-' % color,alpha=alpha,label=label)
		else:
			plt.plot(y*get_sens(data_dir),eyy,'%s-' % color,alpha=alpha)
	#~ fig = plt.figure(data_dir+str(n)+'profile')
	#~ plt.xlim([-20,20])
	plt.grid()
	if xlabel:
		plt.xlabel(r'P coordinate $(mm)$',size=18)
	plt.ylabel(r'Strain field along a generatrix line $\left(mm/mm\right)$',size=18)
	plt.ylim([0.0,ylim])
	if show:
		plt.show()

def plot_tensile_test(data_dir,n=1000,show=False,marker='bo--',lw=1,epsmax=0.):
	file_name = data_dir+get_basename(data_dir)+'-tt.dat'
	step_0,step_f = get_field_params(data_dir)[-2:]
	
	if n < step_f:
		step_f = n
	
	tt = np.loadtxt(file_name)	
	strain = tt[step_0:step_f+1,0]
	stress = tt[step_0:step_f+1,1]*1000/18
	
	if epsmax == 0.:
		epsmax = round(tt[int(get_field_params(data_dir)[-1:][0]),0]*1.05,2)
	
	plt.plot(strain,stress,marker,lw=lw)
	plt.xlabel(r'Strain $(mm/mm)$',size=18)
	plt.ylabel(r'Stress $(MPa)$',size=18)
	
	plt.yticks([0,200,800,1000])
	plt.xticks([0,epsmax])
	plt.grid()
	#~ ax.xaxis.set_label_coords(0.5, -0.025)
	#~ ax.yaxis.set_label_coords(-0.025, 0.5)
	if show:
		plt.show()
	return

def extenso(data_dir,p1_o,p2_o,show=False): #closest node
	def plot_extenso():
		F = get_field(data_dir,'eyy',get_steps(data_dir)[-1])
		fig = plt.figure()
		plt.title('Extensometer : %s->%s' % (p1_o,p2_o))
		plt.plot([p1_o[0],p2_o[0]],[p1_o[1],p2_o[1]],'wo-',lw=2)
		plt.imshow(F,extent=get_extent(data_dir))

		plt.colorbar()
		#~ print str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1])
		#~ print 'extenso--%s,%s-%s,%s' % (str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1]))
		plt.savefig('extenso--%s,%s-%s,%s.png' % (str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1])),bbox_inches='tight',transparent=False,pad_inches=0.1)
		plt.show()
		plt.close()
	if show:
		plot_extenso()
	steps = get_steps(data_dir)

	ext = np.zeros_like(steps)*1.
	
	Lo = np.linalg.norm(p2_o-p1_o)
	u = (p2_o-p1_o)/Lo

	print 'Calculating... '
	for i,step in enumerate(steps):
		F = get_field(data_dir,'eyy',step)

		#~ d1 = np.array([get_value_beta(data_dir,'U',step,p1_o,show=False,verbose=False),get_value_beta(data_dir,'V',step,p1_o,show=False,verbose=False)])
		#~ d2 = np.array([get_value_beta(data_dir,'U',step,p2_o,show=False,verbose=False),get_value_beta(data_dir,'V',step,p2_o,show=False,verbose=False)])
		#~ try:
			#~ d1 = np.array([get_value_beta(data_dir,'U',step,p1_o),get_value_beta(data_dir,'V',step,p1_o)])
			#~ d2 = np.array([get_value_beta(data_dir,'U',step,p2_o),get_value_beta(data_dir,'V',step,p2_o)])
		#~ except:
			#~ d1 = np.array([get_value_beta(data_dir,'u_c',step,p1_o),get_value_beta(data_dir,'v_c',step,p1_o)])
			#~ d2 = np.array([get_value_beta(data_dir,'u_c',step,p2_o),get_value_beta(data_dir,'v_c',step,p2_o)])
		try:
			d1 = np.array([get_value_beta(data_dir,'U',step,p1_o),get_value_beta(data_dir,'V',step,p1_o)])
			d2 = np.array([get_value_beta(data_dir,'U',step,p2_o),get_value_beta(data_dir,'V',step,p2_o)])
		except:
			d1 = np.array([get_value_beta(data_dir,'u_c',step,p1_o),0])
			d2 = np.array([get_value_beta(data_dir,'u_c',step,p2_o),0])

		p1_f = p1_o + d1
		p2_f = p2_o + d2
		
		#~ real extensometer
		Lf = np.linalg.norm(p2_f-p1_f)
		
		#~ extensometer holding initial direction
		#~ Lf = Lo + np.inner(p1_f-p1_o,-u) + np.inner(p2_f-p2_o,u)

		dL = (Lf-Lo)/Lo

		#~ print i, p1_f,p2_f,dL
		

		ext[i] = dL
	print 'Done !'
	
	return ext

def extenso_beta(data_dir,p1,p2):
	def plot_extenso():
		F = get_field(data_dir,'eyy',get_steps(data_dir)[-1])
		fig = plt.figure()
		plt.title('Extensometer : %s->%s' % (p1_o_mm,p2_o_mm))
		plt.plot([p1_o_mm[0],p2_o_mm[0]],[p1_o_mm[1],p2_o_mm[1]],'wo-',lw=2)
		plt.imshow(F,extent=get_extent(data_dir))

		plt.colorbar()
		#~ print str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1])
		#~ print 'extenso--%s,%s-%s,%s' % (str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1]))
		plt.savefig('extenso--%s,%s-%s,%s.png' % (str(p1_o[0]),str(p1_o[1]),str(p2_o[0]),str(p2_o[1])),bbox_inches='tight',transparent=False,pad_inches=0.1)
		plt.show()
		plt.close()
	
	
	steps = get_steps(data_dir)
	step = steps[0]
	
	#~ Here we get the coordinates in pixels of both points
	p1_o = get_close_node(data_dir,p1)
	p2_o = get_close_node(data_dir,p2)

	#~ Here we get the coordinates in mm of both points
	p1_o_mm = np.array([get_value_coord_px(data_dir,'X',step,p1_o),get_value_coord_px(data_dir,'Y',step,p1_o)])
	p2_o_mm = np.array([get_value_coord_px(data_dir,'X',step,p2_o),get_value_coord_px(data_dir,'Y',step,p2_o)])
		
	
	ext = np.zeros(np.shape(steps)[0])
	
	#~ Extensometer initial length
	Lo = np.linalg.norm(p2_o_mm-p1_o_mm)
	#~ Extensometer direction
	u = (p2_o_mm-p1_o_mm)/Lo

	print 'Calculating... ',
	for i,step in enumerate(steps[::]):
		F = get_field(data_dir,'eyy',step)

		d1 = np.array([get_value_coord_px(data_dir,'U',step,p1_o),get_value_coord_px(data_dir,'V',step,p1_o)])
		d2 = np.array([get_value_coord_px(data_dir,'U',step,p2_o),get_value_coord_px(data_dir,'V',step,p2_o)])
		
		p1_f = p1_o_mm + d1
		p2_f = p2_o_mm + d2
		
		#~ real extensometer
		Lf = np.linalg.norm(p2_f-p1_f)
		
		#~ extensometer holding initial direction
		#~ Lf = Lo + np.inner(p1_f-p1_o,-u) + np.inner(p2_f-p2_o,u)

		dL = (Lf-Lo)/Lo

		#~ print i, p1_f,p2_f,dL
		
		ext[i] = dL
		
	#~ plot_extenso()
	print 'Done !'
	return ext

def tensile_test_data_beta(data_dir,ext=25,n=1000,S=18):
	file_name = data_dir+get_basename(data_dir)+'-tt.dat'
	step_0,step_f = get_field_params(data_dir)[-2:]
	
	if n < step_f:
		step_f = n
	
	tt = np.loadtxt(file_name)	
	#~ strain = tt[step_0:step_f,0]
	strain = extenso(data_dir,np.array([0.,-ext/2]),np.array([0.,+ext/2]))
	
	stress = tt[step_0:step_f+1,1]*1000./S
	
	return strain,stress

def plot_tensile_test_beta(data_dir,p1_o,p2_o,n=1000,show=False,marker='bo--',lw=1,epsmax=0.):
	file_name = data_dir+get_basename(data_dir)+'-tt.dat'
	step_0,step_f = get_field_params(data_dir)[-2:]
	
	if n < step_f:
		step_f = n
	
	tt = np.loadtxt(file_name)	
	#~ strain = tt[step_0:step_f,0]
	strain = extenso(data_dir,p1_o,p2_o)
	stress = tt[step_0:step_f+1,1]*1000/18
	
	if epsmax == 0.:
		epsmax = round(tt[int(get_field_params(data_dir)[-1:][0]),0]*1.05,2)
	
	plt.plot(strain,stress,marker,lw=lw)
	plt.xlabel(r'Strain $(mm/mm)$',size=18)
	plt.ylabel(r'Stress $(MPa)$',size=18)
	
	plt.yticks([0,200,800,1000])
	plt.xticks([0,epsmax])
	plt.grid()
	#~ ax.xaxis.set_label_coords(0.5, -0.025)
	#~ ax.yaxis.set_label_coords(-0.025, 0.5)
	if show:
		plt.show()
	return

def plot_space_time_strain(data_dir,T=0.500):
	#~ Getting the steps of the correlation
	steps = get_steps(data_dir)

	#~ Preparing the size in pixels of the time-space plot
	n = np.shape(steps)[0]
	Tx,Ty = get_mesh_size(data_dir)

	#~ Preparing the time-space array 
	data = np.zeros([Ty,n])

	#~ Preparing sizes in time and space scales
	#~ T = 0.500 #s (image samplig every T)
	t0,tf = 0, n*T
	yo,yf = get_extent(data_dir)[:2:]

	#~ Defining the extent for the plot
	extent = t0,tf,yo,yf

	#~ Gathering values
	for i,step in enumerate(steps[::]):	
		eyy=get_field(data_dir,'eyy',step)
		data[::,i] = eyy[::,Tx//2]

	#~ Taking its derivative
	#~ d_data = (data[1::,::]-data[:-1:,::])/T
	
	'''
	fig = plt.figure('d-t')
	plt.imshow(data, cmap='hot', extent=extent, interpolation='nearest',aspect='auto')
	plt.xlabel(r'Time $(s)$',size=20)
	plt.ylabel(r'$y$ coordinate $(mm)$',size=20)
	plt.colorbar()
	plt.show()
	'''
	return data

def get_start_end_steps(data_dir):
	step_0,step_f = get_field_params(data_dir)[-2:]
	return step_0,step_f
'''
def trash ()
#~ def get_mesh(data_dir,basename,steps):
	#~ step=steps[0]
	
	#~ file_name = data_dir +basename +'-' + str(step).zfill(5) +'_0.mat'

	#~ X = np.array((io.loadmat(file_name))['x'])[0,::]
	#~ Y = np.array((io.loadmat(file_name))['y'])[::,0]
	#~ X, Y = np.meshgrid(X, Y)
	#~ return X,Y

	def get_field(weld,file_name,step,field):
		data = np.array([('Ti17-Ti17'		,-0.770, 2,-32,9,-8,0.00),
						 ('Ti6242-Ti6242'	,-0.970,32, -9,9,-8,3.37),
						 ('Ti17-Ti64'		,-1.350,10, -5,5,-4,0.00),
						 ('Ti17-Ti6242'		,-0.975, 3, -2,9,-8,0.00)]
						,dtype=[('weld','S20'),('angle',np.float32),('x0',int),('xf',int),('y0',int),('yf',int),('weld_pos',float)])


		angle		= data[data['weld']==weld]['angle'][0]
		x0			= data[data['weld']==weld]['x0'][0]
		xf			= data[data['weld']==weld]['xf'][0]
		y0			= data[data['weld']==weld]['y0'][0]
		yf			= data[data['weld']==weld]['yf'][0]
		weld_pos	= data[data['weld']==weld]['weld_pos'][0]

		
		Z = np.array((mio.loadmat(file_name))[field])

		rotated=ndimage.rotate(Z, angle, axes=(1, 0), reshape=True)

		rotated_cropped= rotated[x0:xf,y0:yf:]

		if raw_input('Plot?')=='y':
			fig = plt.figure[0]
			plt.subplot(1,3,1)
			plt.title('Raw data')
			plt.imshow(Z,interpolation='nearest')
			plt.axis('off')
			plt.subplot(1,3,2)
			plt.title('Rotated')
			plt.imshow(rotated,interpolation='nearest')
			plt.axis('off')
			plt.subplot(1,3,3)	
	def get_field(weld,file_name,step,field):
		data = np.array([('Ti17-Ti17'		,-0.770, 2,-32,9,-8,0.00),
						 ('Ti6242-Ti6242'	,-0.970,32, -9,9,-8,3.37),
						 ('Ti17-Ti64'		,-1.350,10, -5,5,-4,0.00),
						 ('Ti17-Ti6242'		,-0.975, 3, -2,9,-8,0.00)]
						,dtype=[('weld','S20'),('angle',np.float32),('x0',int),('xf',int),('y0',int),('yf',int),('weld_pos',float)])


		angle		= data[data['weld']==weld]['angle'][0]
		x0			= data[data['weld']==weld]['x0'][0]
		xf			= data[data['weld']==weld]['xf'][0]
		y0			= data[data['weld']==weld]['y0'][0]
		yf			= data[data['weld']==weld]['yf'][0]
		weld_pos	= data[data['weld']==weld]['weld_pos'][0]

		
		Z = np.array((mio.loadmat(file_name))[field])

		rotated=ndimage.rotate(Z, angle, axes=(1, 0), reshape=True)

		rotated_cropped= rotated[x0:xf,y0:yf:]

		if raw_input('Plot?')=='y':
			fig = plt.figure[0]
			plt.subplot(1,3,1)
			plt.title('Raw data')
			plt.imshow(Z,interpolation='nearest')
			plt.axis('off')
			plt.subplot(1,3,2)
			plt.title('Rotated')
			plt.imshow(rotated,interpolation='nearest')
			plt.axis('off')
			plt.subplot(1,3,3)
			plt.title('Rotated and cropped')
			plt.imshow(rotated_cropped,interpolation='nearest')
			plt.axis('off')
			plt.savefig('/home/jgarcia/Documents/these/00-titanium_img/2016-12-06--report_dic/rotation--%s.pdf' % (search+'-'+str(step)+'-'+field), bbox_inches='tight')
			plt.show()
		
		if field=='Y':
			rotated_cropped = rotated_cropped-weld_pos
		return rotated_cropped[::,::-1]
			plt.title('Rotated and cropped')
			plt.imshow(rotated_cropped,interpolation='nearest')
			plt.axis('off')
			plt.savefig('/home/jgarcia/Documents/these/00-titanium_img/2016-12-06--report_dic/rotation--%s.pdf' % (search+'-'+str(step)+'-'+field), bbox_inches='tight')
			plt.show()
		
		if field=='Y':
			rotated_cropped = rotated_cropped-weld_pos
		return rotated_cropped[::,::-1]
def treat_fields(data_dir):
	#~ print data_dir
	print get_field_params(data_dir)
	x,y,angle,pos,sens,dx,dy = get_field_params(data_dir)
	steps = get_steps(data_dir)
	basename = get_basename(data_dir)
	steps = get_steps(data_dir)
	step = steps[0]
	#~ X, Y = get_mesh(data_dir,basename,steps)
	file_name = data_dir +basename +'-' + str(step).zfill(5) +'_0.mat'
	file_out = data_dir+'/treated/'+basename+'-' + str(step).zfill(5) +'_0.mat'
	data = io.loadmat(file_name)

	for field in io.whosmat(file_name):
		#~ plt.subplot(121)
		#~ plt.title('%s %s' % (field[0],np.shape(data[field[0]])))
		#~ plt.imshow(data[field[0]],interpolation='nearest')
		data[field[0]]=ndimage.rotate(data[field[0]], angle, axes=(1, 0), reshape=False)
		data[field[0]]=data[field[0]][x:-x:,y:-y:]
		#~ plt.subplot(122)
		#~ plt.title('%s %s' % (field[0],np.shape(data[field[0]])))
		#~ if field[0]=='X': 
			#~ plt.title('%s %s this is going to be treated' % (field[0],np.shape(data[field[0]])))
			#~ data[field[0]]+=-dx
		
		#~ if field[0]=='Y': 
			#~ plt.title('%s %s this is going to be treated' % (field[0],np.shape(data[field[0]])))
			#~ data[field[0]]+=-dy

		#~ plt.imshow(data[field[0]],interpolation='nearest')
		#~ plt.show()
		io.savemat(file_out,data,appendmat=True)
		#~ io.savemat(file_out,data)
	#~ Z=ndimage.rotate(Z0, angle, axes=(1, 0), reshape=False)
	#~ return X[x:-x:,y:-y:],Y[x:-x:,y:-y:],Z[x:-x:,y:-y:]
	#~ Z[x:-x:,y:-y:]
	return 


def get_fields(data_dir,step): #~ Don't use this yet.
	basename = get_basename(data_dir)
	file_name = data_dir +basename +'-' + str(step).zfill(5) +'_0.mat'
	data = np.array((io.loadmat(file_name)))
	return data


def combined_steps(data_dir1,data_dir2):
	steps1 = get_steps(data_dir1)
	steps2 = get_steps(data_dir2)

	steps12 = np.concatenate([steps1,steps2],axis=0)

	smax = np.max(steps12)+2

	f,b = np.histogram(steps12,range(smax))

	steps = b[f==2]

	return steps
'''

def get_force(data_dir):
	#~ The X-tt.dat file should have the analogic data coming from the load frame. The analogic collected information
	#~ is often the displacement of the head and the force during the tensile test.
	basename = get_basename(data_dir)
	data = np.loadtxt(data_dir+basename+'-tt.dat')
	step_0,stef_f = get_start_end_steps(data_dir)
	return data[step_0:stef_f+1,1]

def fit_circle(a,b,c):
	x = np.array([a[0],b[0],c[0]])
	y = np.array([a[1],b[1],c[1]])
	
	# coordinates of the barycenter
	x_m = np.mean(x)
	y_m = np.mean(y)

	# calculation of the reduced coordinates
	u = x - x_m
	v = y - y_m

	# linear system defining the center (uc, vc) in reduced coordinates:
	#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
	#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
	Suv  = np.sum(u*v)
	Suu  = np.sum(u**2)
	Svv  = np.sum(v**2)
	Suuv = np.sum(u**2 * v)
	Suvv = np.sum(u * v**2)
	Suuu = np.sum(u**3)
	Svvv = np.sum(v**3)

	# Solving the linear system
	A = np.array([ [ Suu, Suv ], [Suv, Svv]])
	B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
	uc, vc = np.linalg.solve(A, B)

	xc_1 = x_m + uc
	yc_1 = y_m + vc

	# Calcul des distances au centre (xc_1, yc_1)
	Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
	R_1      = np.mean(Ri_1)
	residu_1 = np.sum((Ri_1-R_1)**2)

	return R_1,np.array([xc_1,yc_1])

def plot_thumbnail(xc,yc,SS):
	x = np.array([xc-SS/2,xc+SS/2,xc+SS/2,xc-SS/2,xc-SS/2])
	y = np.array([yc-SS/2,yc-SS/2,yc+SS/2,yc+SS/2,yc-SS/2])
	
	plt.plot(x,y,'g-',lw=2)
	plt.plot(xc,yc,'ro',ms=5)
	
	return

def get_thumbnail(img,xc,yc,SS,depth,show=False):
	#~ The subset of size SS is taken around the point of reference.
	#~ if SS % 2 == 1: SS += 1
	d = SS//2
	thumbnail_f = img[yc-d:yc+d+1,xc-d:xc+d+1]

	if depth > 1:
		thumbnail_f	= spline_2D(thumbnail_f,depth)
		#~ zos			= spline_2D(zos,depth)
				
	if show:
		plt.figure(figsize=(20,20))
		plt.subplot(121)
		plt.imshow(img,interpolation='nearest',origin='lower',cmap=cm.gray)
		plt.xlim(0,1600)
		plt.ylim(0,1200)
		plot_thumbnail(xc,yc,SS)
		plt.subplot(122)
		plt.imshow(thumbnail_f,interpolation='nearest',origin='lower',cmap=cm.gray)
		plt.show()
	
	return thumbnail_f

def spline_2D(data,depth=2):
	#~ This function is used for splining botht the thumbnails and
	#~ the zone where a given thumbnail is going searched.
	#~ Degree is the level of resolution magnification that will be used. 
	#~ Standard is two, since this function will increase the splined image
	#~ to the second power of the depth. The splined image will have dimensions
	#~ of (degre x nx) x (degre x ny).
	nx,ny = np.shape(data)
	
	x = np.arange(nx)
	y = np.arange(ny)
	f = interpolate.interp2d(x, y, data, kind='cubic')
	
	
	xx = np.arange(nx,step=1./depth)
	yy = np.arange(ny,step=1./depth)
	data_splined  = f(xx,yy)
	
	return data_splined

def get_zos(img,xc,yc,dmax,SS,depth=1,show=False):
	#~ This function is equivalent to get_thumbnail. It gets the Zone Of Search 
	#~ around a given point (xc,yc) on the g image and a maximal displacement of dmax.
	d=(SS-1)//2+dmax
	zos = img[yc-d:yc+d+1,xc-d:xc+d+1]

	if depth > 1:
		#~ thumbnail_f	= spline_2D(thumbnail_f,depth)
		zos			= spline_2D(zos,depth)

	if show:

		x = np.array([xc-d,xc+d,xc+d,xc-d,xc-d])
		y = np.array([yc-d,yc-d,yc+d,yc+d,yc-d])
	
		

		plt.figure(figsize=(15,15))
		plt.subplot(221)
		vmin= np.min(img)
		vmax= np.max(img)
		plt.imshow(img,interpolation='nearest',origin='lower',cmap=cm.gray,vmin=vmin,vmax=vmax)
		plt.plot(x,y,'g-',lw=2)
		plt.plot(xc,yc,'ro',ms=5)
		plot_thumbnail(xc,yc,SS)
		plt.ylim(0,np.shape(img)[0])
		plt.xlim(0,np.shape(img)[1])

		plt.subplot(222)
		vmin= np.min(img)
		vmax= np.max(img)
		plt.imshow(img,interpolation='nearest',origin='lower',cmap=cm.gray,vmin=vmin,vmax=vmax)
		plt.plot(x,y,'g-',lw=2)
		plt.plot(xc,yc,'ro',ms=5)
		plot_thumbnail(xc,yc,SS)

		plt.xlim(xc-4*d,xc+4*d)
		plt.ylim(yc-4*d,yc+4*d)
				
		
		plt.subplot(223)
		plt.imshow(zos,interpolation='nearest',origin='lower',cmap=cm.gray,vmin=vmin,vmax=vmax)
		for i in range((dmax+SS//2-dmax)*depth,(dmax+SS//2+dmax)*depth+1):
			for j in range((dmax+SS//2-dmax)*depth,(dmax+SS//2+dmax)*depth+1):
				plt.plot(i,j,'kx',mew=.25)
		plt.xlim(-0.5,(SS+dmax*2-0.5)*depth)
		plt.ylim(-0.5,(SS+dmax*2-0.5)*depth)

		plt.subplot(224)
		plt.imshow(zos,interpolation='nearest',origin='lower',cmap=cm.gray,vmin=vmin,vmax=vmax)
		plt.xlim(-0.5,(SS+dmax*2-0.5)*depth)
		plt.ylim(-0.5,(SS+dmax*2-0.5)*depth)
		
		plt.show()
	
	return zos

def get_displ(thumbnail_f,g,xf,yf,dmax,depth=1,name='',show=False):
	#~ This function will find the thumbnail (thumbnail_f) in the deformed 
	#~ image (g) starting at the position (xf,yf) up to a (dmax) px
	#~ distance from the starting point. 
	
	def fill_std(start_x,end_x,start_y,end_y):
		for i in range(start_x,end_x):
			for j in range(start_y,end_y):
				thumbnail_candidat = zos[i:i+SS*depth,j:j+SS*depth]
				std[i,j] = np.sum(np.sum(np.power(thumbnail_candidat-thumbnail_f,2),axis=0),axis=0)
				#~ print(std[i,j])
				
				#~ plt.figure(figsize=(22,22))
				#~ plt.subplot(151)
				#~ plt.imshow(zos,interpolation='nearest',cmap=cm.gray,origin='lower')
				#~ plot_thumbnail((SS*depth)/2+j,(SS*depth)/2+i,SS)
				#~ plt.xlim(0,(SS+2*dmax)*depth)
				#~ plt.ylim(0,(SS+2*dmax)*depth)
				
				#~ plt.subplot(152)
				#~ plt.imshow(thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')
				
				#~ plt.subplot(153)
				#~ plt.title(np.mean(thumbnail_candidat))
				#~ plt.imshow(thumbnail_candidat,interpolation='nearest',cmap=cm.gray,origin='lower')
				
				#~ plt.subplot(154)
				#~ plt.imshow(thumbnail_candidat-thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')
				
				#~ plt.subplot(155)
				#~ plt.imshow(std,interpolation='nearest',origin='lower')
				#~ plt.show()
		return


	SS = np.shape(thumbnail_f)[0]
	zos = get_zos(g,xf,yf,dmax,SS)
	
	
	if depth > 1:
		thumbnail_f	= spline_2D(thumbnail_f,depth)
		zos			= spline_2D(zos,depth)
	
	
	
	std = np.zeros([2*dmax*depth,2*dmax*depth])

	#~ t1 = threading.Thread(target=fill_std, args=(0*dmax*depth,2*dmax*depth,0*dmax*depth,2*dmax*depth))
	t1 = threading.Thread(target=fill_std, args=(0*dmax*depth,1*dmax*depth,0*dmax*depth,1*dmax*depth))
	t2 = threading.Thread(target=fill_std, args=(0*dmax*depth,1*dmax*depth,1*dmax*depth,2*dmax*depth))
	t3 = threading.Thread(target=fill_std, args=(1*dmax*depth,2*dmax*depth,0*dmax*depth,1*dmax*depth))
	t4 = threading.Thread(target=fill_std, args=(1*dmax*depth,2*dmax*depth,1*dmax*depth,2*dmax*depth))
	
	t1.start()
	t2.start()
	t3.start()
	t4.start()
	
	t1.join()
	t2.join()
	t3.join()
	t4.join()	
		
	center = np.unravel_index(std.argmin(), std.shape)
	
	dr = center-np.array([dmax*depth,dmax*depth])
	
	
	thumbnail_candidat = zos[center[0]:center[0]+SS*depth,center[1]:center[1]+SS*depth]
	
	plt.figure(figsize=(22,22))
	plt.subplot(151)
	plt.title('Zone of search')
	plt.imshow(zos,interpolation='nearest',cmap=cm.gray,origin='lower')
	plot_thumbnail((SS*depth)/2+center[1],(SS*depth)/2+center[0],SS*depth)
	plt.xlim(0,(SS+2*dmax)*depth)
	plt.ylim(0,(SS+2*dmax)*depth)
	
	plt.subplot(152)
	plt.title('f thumbnail')
	plt.imshow(thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')
	
	plt.subplot(153)
	plt.title('g thumbnail')
	plt.imshow(thumbnail_candidat,interpolation='nearest',cmap=cm.gray,origin='lower')
	
	plt.subplot(154)
	plt.title('f-g thumbnail')
	plt.imshow(thumbnail_candidat-thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')
	
	plt.subplot(155)
	plt.title('STD thumbnail')
	plt.imshow(std,interpolation='nearest',origin='lower')
	
	plt.savefig('%s.png' % name, bbox_inches='tight',transparent=False,pad_inches=0.1)
	
	if show:
		plt.show()
		
	
	plt.close()
	
	displ = dr*1./depth
	
	return displ

'''
SS = 3
X = np.ones([25,25])*5

def thumb(x,y,SS):
	if SS % 2 == 00: SS+=1
	d = SS//2
	X[y-d:y+d+1,x-d:x+d+1]=2


def zos(x,y,SS,dmax):
	d=SS//2+dmax
	X[y-d:y+d+1,x-d:x+d+1]=1
	xx = np.arange(x-dmax,x+dmax+1)
	yy = np.arange(y-dmax,y+dmax+1)
	
	print dmax*2+1
	print np.shape(xx)
	plt.plot(xx,yy,'wx',mew=2,ms=8)

for dmax in range(8):
	zos(10,10,SS,dmax)
	thumb(10,10,SS)
	plt.imshow(X,interpolation='nearest',origin='lower')
	plt.show()


'''




def correlate(f_path,xf,g_path,xg,SS,dmax,depth,step,show=False,savethumbs=False,savedir=''):
	#~ global g_path
	def plot_correlation(show=False,save=False):
		plt.figure(figsize=(22,22))
		plt.subplot(151)
		plt.title('Zone of search')
		plt.imshow(zos,interpolation='nearest',cmap=cm.gray,origin='lower')
		plot_thumbnail(RF[0],RF[1],SS*depth)
		plt.xlim(0,(SS+2*dmax)*depth)
		plt.ylim(0,(SS+2*dmax)*depth)

		plt.subplot(152)
		plt.title('f thumbnail')
		plt.imshow(thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')

		plt.subplot(153)
		plt.title('g thumbnail')
		plt.imshow(thumbnail_g,interpolation='nearest',cmap=cm.gray,origin='lower')

		plt.subplot(154)
		plt.title('f-g thumbnail')
		plt.imshow(thumbnail_g-thumbnail_f,interpolation='nearest',cmap=cm.gray,origin='lower')

		plt.subplot(155)
		plt.title('STD thumbnail')
		plt.imshow(std[::-1,::-1],interpolation='nearest',origin='lower')
		if dr[1]<>0 and dr[0]<>0:
			plt.arrow(r0[0],r0[1],dr[0],dr[1],color='g',width=.5,length_includes_head=True)
			
		plt.xlim(-.5,np.shape(std)[1]-.5)
		plt.ylim(-.5,np.shape(std)[0]-.5)	
		if savethumbs:
			plt.savefig('%s%s-%.2f-%.f2.png' % (savedir,str(step).zfill(3),xf[0],xf[1]), bbox_inches='tight',transparent=False,pad_inches=0.1)

		if show:
			plt.show()

		plt.close()
		return
	
	def compute_std(dmax,depth):
		std = np.zeros([2*dmax*depth+1,2*dmax*depth+1])
		def fill_std(start_x,end_x,start_y,end_y):
			for i in range(start_x,end_x):
				for j in range(start_y,end_y):
					thumbnail_candidat = zos[i:i+SS*depth,j:j+SS*depth]
					std[i,j] = np.sum(np.sum(np.power(thumbnail_candidat-thumbnail_f,2),axis=0),axis=0)	
					
		#~ fill the std matrix using 4 cores
		#~ t1 = threading.Thread(target=fill_std, args=(0*dmax*depth,2*dmax*depth+1,0*dmax*depth,2*dmax*depth+1))
		t1 = threading.Thread(target=fill_std, args=(0*dmax*depth,1*dmax*depth,0*dmax*depth,1*dmax*depth))
		t2 = threading.Thread(target=fill_std, args=(0*dmax*depth,1*dmax*depth,1*dmax*depth,2*dmax*depth+1))
		t3 = threading.Thread(target=fill_std, args=(1*dmax*depth,2*dmax*depth+1,0*dmax*depth,1*dmax*depth))
		t4 = threading.Thread(target=fill_std, args=(1*dmax*depth,2*dmax*depth+1,1*dmax*depth,2*dmax*depth+1))
		
		t1.start()
		t2.start()
		t3.start()
		t4.start()
		
		t1.join()
		t2.join()
		t3.join()
		t4.join()	
		
		return std	
		
	#~ get the thumbnail to be tracked in the f image
	f = misc.imread(f_path)[::,::,1]
	thumbnail_f = get_thumbnail(f,xf[0],xf[1],SS,depth,show=False)
	
	#~ get zone of search on the g image
	g = misc.imread(g_path)[::,::,1]
	
	zos = get_zos(g,xg[0],xg[1],dmax,SS,depth,show=False)
	
	#~ compute the std of the difference between the reference thumbnail
	#~ and all the thumbnails on the search zone
	std = compute_std(dmax,depth)
	
	#~ find the thumbnail that minimises the std
	rf = np.array(np.unravel_index(std.argmin(), std.shape))[::-1]
	#~ print(rf)
	
	#~ the center of the std array is the reference point
	#~ its center is the barycenter on the f image
	r0 = (np.shape(std)[0]-1)/2.*np.array([1.,1.])
	#~ print(r0)
	
	#~ take the difference of the position of the thumbnail found and the center
	#~ of the std array. This will give the displacement
	dr = (rf-r0)*1.
	#~ print(dr)

	R0 = np.array(np.shape(zos))/2
	RF = R0+dr


	g = misc.imread(g_path)[::,::,1]
	thumbnail_g = get_thumbnail(g,xg[0]+dr[0],xg[1]+dr[1],SS,depth,show=False)
	
	
	plot_correlation(show=show,save=savethumbs)
	
	return dr/depth
		
	

def get_section(name):
	measures = np.loadtxt('/home/jgarcia/Documents/these/12-dic_results/specimens_geometry.dat',usecols=(1,2,3,4))
	names = np.loadtxt('/home/jgarcia/Documents/these/12-dic_results/specimens_geometry.dat',usecols=(0,),dtype='|S8')
	
	mask = (names==name)
	
	t = np.mean(measures[mask,:3:][0])
	W = measures[mask,3][0]
	
	return t*W
