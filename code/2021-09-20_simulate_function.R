## Function to do the simulation:
simulate.data <- function(params,illustration=F){
  ## Convenience functions:
  expit <- function(x) exp(x)/(1+exp(x))
  
  gompertz <- function(drop,slope,x){
    -drop*(exp(-exp(-exp(1)*slope*((x-16)/48))))
  }
  
  ## Put all the param names in local scope 
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  # Vector of market selection indicators:
  s.m=rbinom(N.m,1,p.sel)
  # Number of selected and non-selected markets
  N.s.m = sum(s.m)
  N.n.m = sum(1-s.m)
  
  # Market-level growth params
  Sig.1 = diag(c(sig.1a,sig.1b)) %*% 
    matrix(c(1,rho.1,rho.1,1),2,2) %*% 
    diag(c(sig.1a,sig.1b))
  # Draw market intercepts and slopes
  growth <- MASS::mvrnorm(N.m,c(mu.alpha,mu.beta),Sigma = Sig.1)
  
  # Draw COVID impacts in selected markets
  covid.s <- -1*abs(rnorm(N.s.m,mu.covid,sd.covid))
  recovery.s <- abs(rnorm(N.s.m,mu.recov,sd.recov))
  
  # Draw COVID impacts in non-selected markets
  covid.n <-  -1*abs(rnorm(N.n.m,mu.covid,sd.covid))
  recovery.n <- abs(rnorm(N.n.m,mu.recov,sd.recov))
  
  # For non-selected markets, draw a count of practices 
  N.p.n.m <- round(rexp(N.n.m,1/N.p.n)+1,0) # can't have zero practices
  p.join.n.m <- rep(0,N.n.m)

  # For selected markets, draw a count of practices and a joining probability
  N.p.s.m <- round(rexp(N.s.m,1/N.p.s)+1,0) # can't have zero practices
  p.join.s.m <- abs(rnorm(N.s.m,p.join,.1))
  
  # Collect market-level parameters
  dat.m <-  tibble(id.m=factor(1:N.m),
                   s.m=sort(s.m), # ordered so can concatenate selected and non-selected below
                   int.m=growth[,1],
                   slope.m=growth[,2],
                   N.p.m=c(N.p.n.m,N.p.s.m),
                   p.join.m=c(p.join.n.m,p.join.s.m),
                   drop.covid=c(covid.n,covid.s),
                   recov.covid=c(recovery.n,recovery.s)) %>%
    mutate(selected=factor(s.m,levels=0:1,labels=c('NotSelected','Selected')))
  # Clean up
  rm(s.m,N.n.m,N.p.n.m,p.join.n.m,N.s.m,N.p.s.m,p.join.s.m,Sig.1,growth,covid.s,recovery.s,covid.n,recovery.n)
  
  ## Examine the market-level distributions of slopes
  if (F){
    ggplot(dat.m,aes(x=slope.m)) + geom_density(aes(col=selected))
    dat.m %>% group_by(selected) %>% summarize(count=n())
    dat.m %>% group_by(selected) %>% summarize(mean.slope=mean(slope.m))
  }
  
  # Function to draw practice-level params
  draw.prov.params <- function(df){
    if (df$s.m==1){
      # Var-cov matrix for practices in selected markets
      Sig.2 = diag(c(sig.2a,sig.2b)) %*% 
        matrix(c(1,rho.2,rho.2,1),2,2) %*% 
        diag(c(sig.2a,sig.2b))
      # Draw practice intercepts and slopes
      MASS::mvrnorm(1,c(df$int.m - gamma.1*(1-df$trt.i),
                        df$slope.m - gamma.2*(1-df$trt.i)),Sigma = Sig.2)
    } else {
      # Var-cov matrix for practices in non-selected markets
      Sig.2 = diag(c(sig.2a,sig.2b)) %*% 
        matrix(c(1,rho.2,rho.2,1),2,2) %*% 
        diag(c(sig.2a,sig.2b))
      # Drawn practice intercepts and slopes
      MASS::mvrnorm(1,c(df$int.m - theta.1,
                        df$slope.m - theta.2),Sigma = Sig.2)
    }
  }
  
  # Assemble practice-level data and params
  dat.p <- as_tibble(data.frame(id.m=factor(rep(dat.m$id.m,dat.m$N.p.m)),
                                id.i=factor(unlist(sapply(dat.m$N.p.m,seq,from=1))),
                                p.join.i=rep(dat.m$p.join.m,dat.m$N.p.m))) %>% 
    left_join(dat.m,by='id.m') %>%  rowwise() %>%
    # Simulate joining and matching
    mutate(trt.i = ifelse(s.m==1,rbinom(1,1,p.join.i),0),
           match.i = ifelse(trt.i==0,rbinom(1,1,p.match),0)) %>% ungroup() %>%
    group_by(id.m,id.i) %>% nest() %>%
    mutate(params=map(data,draw.prov.params)) %>%
    unnest(cols=c(data,params)) %>% mutate(name=c('int.i','slope.i')) %>% 
    pivot_wider(names_from=name,values_from=params) %>% ungroup() %>%
    mutate(class=factor(trt.i,levels=c(0,1),labels=c('Ctrl','PCF')))
  
  ## Check the "lumpiness"
  if (F) {
    temp <- dat.p %>% group_by(selected,id.m,trt.i) %>% summarize(count=n()) %>%
      ungroup() 
    ggplot(filter(temp,selected=="NotSelected"),aes(x=count)) + geom_histogram() +
      ggtitle("Non-Selected Markets: Count of Practices")
    ggplot(filter(temp,selected=="Selected") %>% 
             pivot_wider(names_from=trt.i,values_from=count) %>%
             mutate(ctrl=ifelse(is.na(`0`),0,`0`),
                    trt=ifelse(is.na(`1`),0,`1`)),
           aes(x=ctrl,y=trt)) + geom_point() + xlab('Non-Joiners') + ylab('Joiners') +
      ggtitle("Selected Markets: Count of Practices") + geom_abline(slope=.16)
    
    dat.p %>% group_by(s.m,trt.i,match.i) %>% summarize(count=n())
  }
  
  # Make a timeseries for each provider: 24 months of pre and 12*5=60 months of post
  full.dat <- tibble(month=rep(1:84,sum(dat.m$N.p.m)),
                     id.m=rep(dat.p$id.m,each=84),
                     id.i=rep(dat.p$id.i,each=84)) %>%
    mutate(id.m=factor(id.m),id.i=factor(id.i))
  
  dat.p.t.prelim <- left_join(dat.p,full.dat,by=c("id.m","id.i")) %>%
    mutate(spend.nocovid = int.i + slope.i*month)
  
  # Isolate Feb 2020 as reference year for percent change COVID impacts:
  dat.p.ref <- filter(dat.p.t.prelim,month==14) %>% rename(ref.spend=spend.nocovid) %>%
    select(id.m,id.i,ref.spend)
  # Generate spending in three periods: before April 2020, April 2020, After April 2020
  dat.p.time <- left_join(dat.p.t.prelim,dat.p.ref) %>%
    mutate(spend.covid = ifelse(month<16,spend.nocovid,
                                ifelse(month==16,spend.nocovid+(ref.spend/100*drop.covid),
                                       ifelse(month>16,spend.nocovid+(ref.spend/100)*(drop.covid+gompertz(drop.covid,recov.covid,month)),
                                              NA))))
  
  ## Plot some to check
  if (F){
    ggplot(filter(dat.p.time,id.m%in%c(1:3,304:306)),aes(x=month,y=spend.covid,group=id.i)) + 
      geom_line(aes(col=trt.i)) + facet_wrap(~id.m)
  }
  ## Clean up
  rm(full.dat,dat.p.t.prelim,dat.p.ref)

  if (illustration) {
    return(dat.p.time)
  } else {
    # Provider-level data on analytical sample (i.e., treated and matched controls)
    prov.prelim <- dat.p.time %>% 
      # Take only treated and matched control practices
      filter((trt.i==1)|(match.i==1)) %>%
      mutate(period=ifelse(month<=24,'Pre',
                           paste('Post',((month-25)%/%12)+1,sep="_")))
    # Summarizing the counts of providers in each group
    prov.counts <- filter(prov.prelim,month==1) %>% 
      group_by(selected,class,period) %>% summarize(prov.count=n()) %>%
      select(selected,class,prov.count)
      
    # Summarizing the annual post vs pre changes for covid and no covid scenarios  
    prov.changes <- prov.prelim %>% 
      group_by(selected,class,period) %>%
      summarize(mean.nocovid=mean(spend.nocovid),mean.covid=mean(spend.covid)) %>%
      pivot_longer(mean.nocovid:mean.covid,names_to='covid') %>%
      unite('periodcovid',period:covid) %>%
      pivot_wider(names_from=periodcovid,values_from=value) %>%
      mutate(covid_1=Post_1_mean.covid-Pre_mean.covid,
             nocovid_1=Post_1_mean.nocovid-Pre_mean.nocovid,
             covid_2=Post_2_mean.covid-Pre_mean.covid,
             nocovid_2=Post_2_mean.nocovid-Pre_mean.nocovid,
             covid_3=Post_3_mean.covid-Pre_mean.covid,
             nocovid_3=Post_3_mean.nocovid-Pre_mean.nocovid,
             covid_4=Post_4_mean.covid-Pre_mean.covid,
             nocovid_4=Post_4_mean.nocovid-Pre_mean.nocovid,
             covid_5=Post_5_mean.covid-Pre_mean.covid,
             nocovid_5=Post_5_mean.nocovid-Pre_mean.nocovid) %>%
      select(class,selected,covid_1:nocovid_5) %>% 
      left_join(prov.counts) %>%
      pivot_longer(covid_1:prov.count,names_to = 'covid',values_to='delta') %>%
      unite(col='scenario',class,selected,covid) %>% 
      pivot_wider(names_from=scenario,values_from=delta)
    return(unlist(prov.changes))
  }
}

