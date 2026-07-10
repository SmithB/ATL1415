#! /usr/bin/env python

from scipy.interpolate import interpn
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc

class dzPicker(object):
    def __init__(self, img_data=None, img_args=None, fig=None,
                 field='dz', time_field='t', time_dim=2,
                 handles=None, file_dict=None, dz_dict=None, 
                 file_args=None, W=2.e3):
        if fig is None:
            fig, handles = self.__init_ui__(img_data, img_args)
        
        self.handles=handles
        self.dz_dict=dz_dict
        self.time_field=time_field
        self.field=field
        self.messages=[[]]
        self.last_pt=[[]]
        self.file_dict=file_dict
        if file_args is None:
            self.file_args={}
        else:
            self.file_args=file_args
        self.dz_dict=dz_dict
        self.W=W
        self.cid = fig.canvas.mpl_connect('button_press_event', self)
    
    def __init_ui__(img_data, img_args=None):
        fig=plt.figure()
        hax=fig.subplots(1,2)
        handles={'map_ax':hax[0], 'plot_ax':hax[1]}
        img_data.show(ax=handles['map_ax'], **img_args)
        return(fig, handles)
    
    def __call__(self, event):
        try:
            xy0=(event.xdata, event.ydata)
            tx = 'xy =[%f,%f]' % xy0
            self.handles['plot_ax'].set_title(tx)
            if self.dz_dict is not None:
                dz_dict=self.dz_dict
            elif self.file_dict is not None:
                dz_dict={}
                for key, file in self.file_dict.items():
                    pad=np.array([-0.5, 0.5])*self.W
                    dz_dict[key]=pc.grid.data().from_h5(file, bounds=[xy0[0]+pad, xy0[1]+pad], **self.file_args)
            for key, dz0 in dz_dict.items():
                tt=getattr(dz0, self.time_field)
                self.last_pt += [[key]]
                if self.time_dim==2:
                    zz=interpn((dz0.y, dz0.x, dz0.t), getattr(dz0, self.field), 
                               (event.ydata*np.ones_like(tt), event.xdata*np.ones_like(tt), tt))
                else:
                    zz=interpn((t, dz0.y, dz0.x), getattr(dz0, self.field), 
                               (tt, event.ydata*np.ones_like(tt), event.xdata*np.ones_like(tt))
                self.handles['plot_ax'].plot(tt, zz, label=tx+' '+str(key))
            y_vals=np.r_[[item._y.ravel() for item in self.handles['plot_ax'].lines]].ravel()
            self.handles['plot_ax'].set_ylim([np.nanmin(y_vals), np.nanmax(y_vals)])
        except Exception as e:
            self.messages += [e]
            plt.gca().set_title('ERROR')
        self.handles['plot_ax'].figure.canvas.draw()
    
    def clear_lines(self):
        lines=list(self.handles['plot_ax'].lines)
        for line_no in range(len(list(self.handles['plot_ax'].lines))):
            self.handles['plot_ax'].lines.pop(0)
        self.handles['plot_ax'].figure.canvas.draw()

