plot_native_ternary <- function(stats,out) {
  stats <- stats[match(c("PCOS-CAM B","ScType","SingleR"),stats$method),]
  stopifnot(!anyNA(stats$method), all(stats$fine_clusters<=stats$major_match))
  p <- cbind(other=stats$reviewed-stats$major_match,
             broad=stats$major_match-stats$fine_clusters,
             finer=stats$fine_clusters)/stats$reviewed
  stopifnot(all(abs(rowSums(p)-1)<1e-10))
  H <- sqrt(3)/2
  xy <- function(a,b,c) c(b+c/2,H*c)
  seg <- function(a,b,...) segments(a[1],a[2],b[1],b[2],...)
  gold <- "#F1BD59"; ink <- "#181818"
  draw <- function() {
    par(mar=rep(0,4),family="sans",fg=ink,col=ink,xaxs="i",yaxs="i")
    plot.new();plot.window(c(-.19,1.91),c(-.22,1.15),asp=1)
    text(-.165,1.09,"Ovarian cell subtype annotation",adj=0,font=2,cex=2.25)
    text(-.163,1.017,"GSE268919   /   10 reviewed clusters   /   Each tool uses its own reference",adj=0,cex=1.16)
    polygon(c(0,1,.5),c(0,0,H),col="#ECECEB",border=NA)
    for(t in seq(.1,.9,.1)){
      seg(xy(t,1-t,0),xy(t,0,1-t),col="white",lwd=2)
      seg(xy(1-t,t,0),xy(0,t,1-t),col="white",lwd=2)
      seg(xy(1-t,0,t),xy(0,1-t,t),col="white",lwd=2)
    }
    polygon(c(0,1,.5),c(0,0,H),border=ink,lwd=2.4)
    polygon(c(.25,.75,.5),c(H/2,H/2,0),border="#626262",lwd=1.3)
    for(t in seq(.2,.8,.2)){
      text(t,-.030,round(100*t),cex=.9,font=2)
      text(1-t/2+.041,H*t,round(100*t),cex=.9,font=2)
      text((1-t)/2-.041,H*(1-t),round(100*t),cex=.9,font=2)
    }
    text(.5,H+.055,"Family match + finer call",font=2,cex=1.28)
    text(0,-.074,"No full family\nmatch",cex=1.13,font=2)
    text(1,-.074,"Family match\nonly",cex=1.13,font=2)
    shapes <- c(21,22,24)
    for(i in seq_len(3)) {
      v <- xy(p[i,1],p[i,2],p[i,3])
      points(v[1],v[2],pch=shapes[i],bg=gold,col=ink,lwd=1.8,cex=2.65)
    }
    x0 <- 1.13
    text(x0,.79,"Distinct finer labels",adj=0,font=2,cex=1.6)
    text(x0,.718,"Within matching families",adj=0,cex=1.12)
    text(1.61,.594,"Called",font=2,cex=1.05)
    text(1.81,.594,"Supported*",font=2,cex=1.05)
    segments(x0,.55,1.90,.55,lwd=1.5)
    for(i in seq_len(3)){
      y <- .466-(i-1)*.116
      points(x0+.024,y,pch=shapes[i],bg=gold,col=ink,lwd=1.5,cex=1.85)
      text(x0+.067,y,stats$method[i],adj=0,font=2,cex=1.22)
      text(1.61,y,stats$distinct_fine[i],font=2,cex=1.45)
      text(1.81,y,stats$supported[i],font=2,cex=1.45)
    }
    segments(x0,.161,1.90,.161,lwd=1.5)
    text(x0,.095,"* Local reference supports the specific identity.",adj=0,cex=.85)
    text(-.16,-.173,"Each point represents one method. The three cluster percentages sum to 100%.",adj=0,cex=.96)
  }
  pdf(file.path(out,"native_subtype_comparison.pdf"),width=13,height=8.5,useDingbats=FALSE)
  draw();dev.off()
  png(file.path(out,"native_subtype_comparison.png"),width=13,height=8.5,units="in",res=300,bg="white",type="cairo")
  draw();dev.off()
}
