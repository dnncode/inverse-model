#Scaled form ADMM
# 1). with close-form solutions
# 2). with high frequencies terms are forced to zero
# 3). Lasso + Fused Tikhonov

#y = y.ext; x = x.ext; s=s.ext; r=1; lambda=10
fused.L2 = function(y, x, r, s, lambda1, lambda2, q.high.case){
  
  xx = x[,c(q.high.case, q.high.case+ncol(x)/2)]
  ss = s[,c(q.high.case, q.high.case+ncol(s)/2)]
  w = r 
  
  # initialization
  set.seed(2)
  z1.old = z2.old = u.old = rnorm(2*length(q.high.case),0,0.1)
  
  
  # iteration
  judge = 0
  k = 0
  while (judge<1){
    k = k + 1
    
    # z1 loop -------------
    tmp1 = t(xx)%*%xx + r/2 * diag(ncol(xx))
    tmp2 = t(xx) %*% y + (z2.old - u.old)*r
    z1.new = chol2inv(chol(tmp1)) %*% tmp2
    
    # z2 loop -------------  
    judge2 = 0
    k2 = 0
    z1.old.inner = z2.old.inner = u.old.inner = rnorm(length(z2.old),0,0.1)
    while (judge2<1){
      k2 = k2 + 1
      
      tmp1 = t(ss)%*%ss + (r+w)/2 * diag(ncol(ss))
      tmp2 = z1.new + u.old + z2.old.inner*w - w*u.old.inner
      z1.new.inner = r/2/lambda2 * chol2inv(chol(tmp1)) %*% tmp2
      
      tmp = z1.new.inner + u.old.inner
      z2.new.inner = array(0,dim=c(length(z2.old.inner),1))
      case = which(tmp > lambda1/w) 
      z2.new.inner[case] = tmp[case] - lambda1/w
      case = which(tmp < -lambda1/w) 
      z2.new.inner[case] = tmp[case] + lambda1/w

      u.new.inner = u.old.inner + z1.new.inner - z2.new.inner
      
      # objective function: inner loop
      part1 =  lambda2 * sum( z1.old.inner * ( (t(ss)%*%ss + r/2)%*%z1.old.inner) )
      part2 = r * sum( z1.old.inner *  (u.old + z1.new)  )
      part3 = lambda1 * sum( abs(z2.old.inner) )
      obj.old = part1 + part2 + part3
      
      part1 =  lambda2 * sum( z1.new.inner * ( (t(ss)%*%ss + r/2)%*%z1.new.inner) )
      part2 = r * sum( z1.new.inner * (u.old + z1.new)  )
      part3 = lambda1 * sum( abs(z2.new.inner) )
      obj.new = part1 + part2 + part3
      
      if ( (abs(obj.new-obj.old)<0.1)+(abs(obj.new-obj.old)/abs(obj.old)<0.0001) >= 1){
        judge2 = 1
      }
      
      z1.old.inner = z1.new.inner
      z2.old.inner = z2.new.inner
      u.old.inner = u.new.inner
      
    }
    z2.new = z2.new.inner
    
    # u loop
    u.new = u.old +  (z1.new - z2.new)
    
    # objective function
    part1 =  sum( (y - xx%*%z1.old)^2 )/2
    part2 =  lambda2 * sum( ( ss %*% z2.old )^2 )
    part3 =  lambda1 * sum( abs(z2.old) )
    obj.old = part1 + part2 + part3
    
    part1 =  sum( (y - xx%*%z1.new)^2 )/2
    part2 =  lambda2 * sum( ( ss %*% z2.new )^2 )
    part3 =  lambda1 * sum( abs(z2.new) )
    obj.new = part1 + part2 + part3
    
    if ( (abs(obj.new-obj.old)<0.1)+(abs(obj.new-obj.old)/abs(obj.old)<0.0001) >= 1){
      judge = 1
    }
    
    z1.old = z1.new
    z2.old = z2.new
    u.old = u.new
    print(c(k,obj.old,obj.new))
    
  }
  
  output.list = list()
  output.list[[1]] = z1.new
  output.list[[2]] = z2.new
  return(output.list)
}

