
############################# LOAD PACKAGES ###############################
require(fda)

############################# LOAD FUNCTIONS ##############################

trapzc = function(t_step,y)
{ # compute the integral of y with step "t_step", with trapezoid rule
  return(t_step*(0.5*y[1]+sum(y[2:(length(y)-1)]) +0.5*y[length(y)]))
}

clr = function(density, z, z_step)
{ # transform a density to a clr
  return(log(density)-trapzc(z_step,log(density))/(max(z)-min(z)))
}

clr2density <- function(clr, z, z_step)
{ # back-transform a clr to a density
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

d.A <- function(z,z_step,y1,y2)
{ #Compute the distanze between two densities y1 and y2 defined on the discretized equispaced interval z with step z_step
  x1=clr(z, z_step, y1)
  x2=clr(z, z_step, y2)
  return(sqrt(trapzc(z_step, (x1-x2)^2)))
}

bpol.01 = function(m,x,k,a=0,b=1)
{ # Bernstein polynomials
  y=sc.01(x,a,b)
  return(choose(m,k)*y^k*(1-y)^(m-k))
  #return(m*k*y^k*(1-y)^(m-k))
}

Fd.01 = function(Fn, x, a=0, b=1)
{ # Evaluates the bernstein distribution function given the empirical distribution function, with support in [a,b]
  n=length(knots(Fn))
  m=n
  y=scInv.01(x,a,b)
  ret=rep(0,length(y))
  for(k in 0:m)
  {
    bp=choose(m,k)*y^k*(1-y)^(m-k)
    ret=apply(rbind(ret,Fn(sc.01(k/m,a,b))*bp),2,sum)    
  }
  return(ret)
}

fd.01 = function(Fn, x, a=0, b=1)
{  # Evaluates the bernstein distribution function given the empirical distribution function, with support in [a,b]
  n=length(knots(Fn))
  m=n
  y=scInv.01(x,a,b)
  ret=rep(0,length(y))
  for(k in 0:(m-1))
  {
    bp=choose(m-1,k)*(y^k)*(1-y)^(m-1-k)
    #bp=(m-1)*k*(y^k)*(1-y)^(m-1-k)
    ret=apply(rbind(ret,(Fn(sc.01((k+1)/m,a,b))-Fn(sc.01(k/m,a,b)))*bp),2,sum)
  }
  ret=m*ret*dscInv.01(x,a,b)
  return(ret)
}


sc.01=function(y,a,b)
{ # map y from (a,b) to (0,1)
  if(a==-Inf & b==Inf)
    return(tan(pi*y-pi/2))
  if(a==-Inf & b!=Inf)
    return(((1+b)*y-1)/y)  
  if(a!=Inf & b!=Inf)
    return((b-a)*y+a)
  print('Cannot compute the transformation, support incorrect')
  return(-1)
}


scInv.01=function(x,a,b)
{ # map y from (0,1) to (a,b)
  if(a==-Inf & b==Inf)
    return(.5+atan(x)/pi)
  if(a==-Inf & b!=Inf)
    #return((b-x)/(1+b-x))
    return(1/(1+b-x))
  if(a!=Inf & b!=Inf)
    return((x-a)/(b-a))
  print('Cannot compute the transformation, support incorrect')
  return(-1)
}

dscInv.01=function(x,a,b)
{ # Derivation of change support transformation
  if(a==-Inf & b==Inf)
    return(1/(pi*(1+x^2)))
  if(a==-Inf & b!=Inf)
    #return((b-x)/(1+b-x))
    return(1/(1+b-x)^2)
  if(a!=Inf & b!=Inf)
    return(1/(b-a))
  print('Cannot compute the transformation, support incorrect')
  return(-1)
}
