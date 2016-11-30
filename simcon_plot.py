import argparse as arp
import numpy as np
import re,os,sys

p=arp.ArgumentParser(prog='asciiPlot',description='Script to plot data from csv type files, it accepts any separator and extension')
p.add_argument('--version',action='version',version='%(prog)s 2.4.3')
p.add_argument('files',metavar='FILES',help='Name of the csv type file/s with the extension',nargs='+')
p.add_argument('-c','--columns',help='Columns to be plotted, the first one will be used as x values. The separation between column numbers for the same file should be a comma and between different files a space. To indicate that a column has error bars it can be indicated entering the column as <col>e<errcol>, i. e. the way to indicate that the second column should be the X values, the third the Y values and the 5th is the error of the 3rd, would be ''-co 2,3e5''',nargs='+')
p.add_argument('-mp',help='Do not plot with matplotlib',action='store_false',default=True)
p.add_argument('-eps',help='Write plot as an Encapsulated PostScript. The name can be specified',type=str,nargs='?',default='noeps',const='sieps')
p.add_argument('-xg',help='Plot with xmgrace',action='store_true',default=False)
p.add_argument('-sh','--shade',help='Shade the area between errorbars instead of showing errorbars. Transparency can be specified (float between 0-1). Note: all error areas will have the same transparency',type=float,nargs='?',default=(-1),const=0.5)
p.add_argument('-l','--legend',help='Labels for the legend. Default text is name of file without extension',nargs='+')
p.add_argument('-t','--title',help='Title of the plot',type=str,default='')
p.add_argument('-lt','--legtit',help='Title for the legend',type=str,default='')
p.add_argument('-xl','--xlabel',help='Label of the x axis',type=str,default='')
p.add_argument('-yl','--ylabel',help='Label of the y axis',type=str,default='')
p.add_argument('-x','--xaxis',help='Set the xaxis limits',nargs=2,type=float)
p.add_argument('-y','--yaxis',help='Set the yaxis limits',nargs=2,type=float)
p.add_argument('-sc','--scale',help='Choose the scale between: liner, semilogy, semilogx or logxy',type=str,choices=['lin','logy','logx','logxy'])
p.add_argument('-s','--separator',help='Separator character/s or regular expression to match the separator',type=str)
p.add_argument('-co','--colors',help='Colors to be used in the plot, they must be valid matplotlib colors',nargs='+')
args=p.parse_args()
l1=len(args.files)
colo=['blue','black','red','green','grey','navy','cyan','slategrey']*3
colus=['0,1']*l1

if args.columns:
    ll=0
    for col in args.columns:
        colus[ll]=col
        ll+=1
if args.colors:
    ll=0
    for col in args.colors:
        colo[ll]=col
        ll+=1

if (args.mp or args.eps!='noeps'):
    import matplotlib.font_manager as fnt
    import matplotlib.pyplot as plt
    xlab=args.xlabel
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

colcount=0
legcount=0
fpat=re.compile(r'(?P<nom>[^/\.]+)\.')
errpat=re.compile(r'([0-9]+)e([0-9]+)')
for tiq,doc in enumerate(args.files):
    docols=[ccc for ccc in re.split(',',colus[tiq])]
    coltup=()
    maxcol=docols[:]
    for ll,ccc in enumerate(docols):
        docols[ll]=ccc
        errorplot=errpat.match(ccc)
        if errorplot:
            maxcol[ll]=max([int(ee) for ee in errorplot.groups(1)])
            coltup+=(errorplot.groups(1)[0], errorplot.groups(1)[1])
        else:
            maxcol[ll]=int(ccc)
            coltup+=(int(ccc), )

    numplots=len(docols)-1
    if args.legend:
        leg=args.legend[legcount:legcount+numplots]
        legcount+=numplots
    else:
        seed=fpat.search(doc)
        if seed:
            leg=[str(seed.group('nom'))]*numplots
        else:
            leg=['']*numplots
    if args.separator:
        data=np.loadtxt(doc,dtype=float,delimiter=args.separator)
    else:
        data=np.loadtxt(doc,dtype=float)        
    x=data[:,int(docols[0])]
    if (args.mp or args.eps!='noeps' or args.xg):
        for i in xrange(numplots):
            cc=docols[i+1]
            errorplot=errpat.match(cc)
            if errorplot:
                y=data[:,int(errorplot.groups(1)[0])]
                eror=data[:,int(errorplot.groups(1)[1])]
                if (args.mp or args.eps!='noeps'):
                    if args.shade==-1:
                        graf.errorbar(x,y,yerr=eror,color=colo[colcount],label=leg[i])
                    else:
                        graf.plot(x, y, color=colo[colcount],label=leg[i])
                        graf.fill_between(x,y-eror,y+eror,alpha=args.shade,facecolor=colo[colcount],edgecolor='k',linewidth=0)
                if args.xg:
                    df=xg.DataXYDY(x=x,y=y,dy=eror,symbol=s1,line=l1)
            else:
                if (args.mp or args.eps!='noeps'):
                    graf.plot(x,data[:,int(cc)],color=colo[colcount],label=leg[i])
                if args.xg:
                    df=xg.Data(x=x,y=data[:,int(cc)],symbol=s1,line=l1)
            colcount+=1

if (args.mp or args.eps!='noeps'):
    graf.set_xlabel(xlab)
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
