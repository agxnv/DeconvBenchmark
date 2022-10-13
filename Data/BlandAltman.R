BA.plot <- function(df, title = NULL, xlab = FALSE, plot=FALSE, scale = NULL, graph = NULL){
  df.summary <- ddply(df, .(Method), summarise, mean = mean(difs, na.rm=T), sd = sd(difs, na.rm=T), 
                      min= min(betahat,na.rm=T), max = max(betahat,na.rm=T), cor = cor(betahat, betasim))
  qq <-c(0.1,NA, 0.8*max(df$betasim))
  qq[2] <- qq[1]+ (qq[3]-qq[1])/2
  qq <- round(qq,2)
  if(graph == "LM22"){
    p <- ggplot(df, aes(betasim, difs)) + 
      geom_point(size =0.5) + 
      geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
      geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
      geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
      geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
      geom_smooth() +  scale_x_continuous(breaks = qq)+ 
      theme_bw() + theme(axis.title.x = element_text(size=22,colour = "white"),
                         axis.title.y = element_text(size=22),
                         axis.text.x.bottom = element_text(size=22),
                         axis.text.y.left = element_text(size=22),
                         plot.title = element_text(hjust = 0.5,size = 30),
                         axis.title = element_text(size = 30),
                         strip.text = element_text(size = 26),
                         panel.grid = element_blank())+
      ylim(-.45, 0.6)
  }else{
    p <- ggplot(df, aes(betasim, difs)) + 
      geom_point(size =0.5) + 
      geom_hline(yintercept = 0, colour = "blue",linetype = "dotted") +
      geom_hline(data=df.summary,aes(yintercept=round(mean,3)), color = "red", linetype = "solid") +
      geom_hline(data=df.summary,aes(yintercept=round(mean+2*sd,3)), color = "red", linetype = "dashed") + 
      geom_hline(data=df.summary,aes(yintercept=round(mean-2*sd,3)), color = "red", linetype = "dashed") +
      geom_smooth() +  scale_x_continuous(breaks = qq)+ 
      theme_bw() + theme(axis.title.x = element_text(size=22,colour = "white"),
                         axis.title.y = element_text(size=0),
                         axis.text.x.bottom = element_text(size=22),
                         axis.text.y.left = element_text(size=22),
                         plot.title = element_text(hjust = 0.5,size = 30),
                         axis.title = element_text(size = 30),
                         strip.text = element_text(size = 26),
                         panel.grid = element_blank())+
      ylim(-.45, 0.6)}
  p + ggtitle(title)+
    if(!is.null(scale)) p <- p + scale
  if(!is.null(title)) p <- p + ggtitle(title)
  if(xlab) {
    p <- p + labs( x = "Simulated coefficients.",y = "Error")
  }else p <- p + labs( x = "",y = "Error")
  p <-  p +    facet_wrap(~Method, nrow = 1)
  if(plot) print(p)
  return(p)
}