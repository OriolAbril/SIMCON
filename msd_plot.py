import argparse as arp
import numpy as np
import re,os,sys

p=arp.ArgumentParser(prog='gdrPlot',description='Script to plot data from .dat files obtained with water.f')
p.add_argument('--version',action='version',version='%(prog)s 1.0')
p.add_argument('files',metavar='FILES',help='Name of the csv type file/s with the extension',nargs='+')
p.add_argument('-mp',help='Do not plot with matplotlib',action='store_false',default=True)
p.add_argument('-eps',help='Write plot as an Encapsulated PostScript. The name can be specified',type=str,nargs='?',default='noeps',const='sieps')
p.add_argument('-xg',help='Plot with xmgrace',action='store_true',default=False)
p.add_argument('-l','--legend',help='Labels for the legend. Default text is name of file without extension',nargs='+')
p.add_argument('-t','--title',help='Title of the plot',type=str,default='')
p.add_argument('-lt','--legtit',help='Title for the legend',type=str,default='')
p.add_argument('-xl','--xlabel',help='Label of the x axis',type=str,default='')
p.add_argument('-yl','--ylabel',help='Label of the y axis',type=str,default='')
p.add_argument('-x','--xaxis',help='Set the xaxis limits',nargs=2,type=float)
p.add_argument('-y','--yaxis',help='Set the yaxis limits',nargs=2,type=float)
p.add_argument('-sc','--scale',help='Choose the scale between: liner, semilogy, semilogx or logxy',type=str,choices=['lin','logy','logx','logxy'])
p.add_argument('-co','--colors',help='Colors to be used in the plot, they must be valid matplotlib colors',nargs='+')
args=p.parse_args()
l1=len(args.files)
colo=['blue','black','red','green','grey','navy','cyan','slategrey']*3
dt=5e-4 #picosegons

if args.colors:
    for index,col in enumerate(args.colors):
        colo[index]=col

if (args.mp or args.eps!='noeps'):
    import matplotlib.font_manager as fnt
    import matplotlib.pyplot as plt
    fig=plt.figure(1)
    graf=fig.add_axes([0.13,0.1, 0.8, 0.8])

if args.xg:
    scriptpath=os.path.dirname(os.path.realpath(__file__))
    sys.path.append(scriptpath+'/graceplot')
    import GracePlot as xg
    pgr=xg.GracePlot()
    pg=pgr[0]
    pg.title(args.title)
    s1=xg.Symbol(symbol=0,fillcolor=0)
    l1=xg.Line(type=1,linewidth=1)
    xlabel=xg.Label(args.xlabel)
    ylabel=xg.Label(args.ylabel)
    pg.xaxis(label=xlabel)
    pg.yaxis(label=ylabel)
    if (args.scale=='logx' or args.scale=='logxy'):
        pg.xaxis(scale='logarithmic')
    if (args.scale=='logy' or args.scale=='logxy'):
        pg.yaxis(scale='logarithmic')

fpat=re.compile(r'(?P<nom>[^/\.]+)\.')
for tiq,doc in enumerate(args.files):
    if args.legend:
        leg=args.legend[tiq]
    else:
        seed=fpat.search(doc)
        if seed:
            leg=str(seed.group('nom'))
        else:
            leg=''
    f=open(doc,'r')
    line=f.readlines()[1]
    f.close()
    msd=np.array(re.split(r'\s+',line)[1:-1],dtype=float)
    x=np.arange(0,len(msd))*dt
    if (args.mp or args.eps!='noeps' or args.xg):
        if (args.mp or args.eps!='noeps'):
            graf.plot(x,msd,color=colo[tiq],label=leg)
        if args.xg:
            df=xg.Data(x=x,y=msd,symbol=s1,line=l1)

if (args.mp or args.eps!='noeps'):
    graf.set_xlabel(args.xlabel)
    graf.set_ylabel(args.ylabel)
    if args.yaxis: 
        graf.set_ylim(args.yaxis)
    if args.xaxis:
        graf.set_xlim(args.xaxis)
    if (args.scale=='logx' or args.scale=='logxy'):
        graf.set_xscale('log')
    if (args.scale=='logy' or args.scale=='logxy'):
        graf.set_yscale('log')
    fig.suptitle(args.title)
    graf.grid(True)
    graf.legend(loc='best', prop=fnt.FontProperties(size='medium'),title=args.legtit)
if args.xg:
    pg.plot(df)
    if args.xaxis:
        pg.xaxis(xmin=args.xaxis[0], xmax=args.xaxis[1])
    if args.yaxis:
        pg.yaxis(ymin=args.yaxis[0], ymax=args.yaxis[1])
if args.eps!='noeps':
    if args.eps!='sieps':
        seed=fpat.search(args.eps)
        if seed:
            figname=seed.group('nom')+'.eps'
        else:
            figname=args.eps+'.eps'
    else:
        seed=fpat.search(args.files[0])
        if seed:
            figname=seed.group('nom')+'.eps'
        else:
            figname=args.files[0]+'.eps'
    fig.savefig(figname, format='eps', dpi=1000)
if args.mp:    
    plt.show()
