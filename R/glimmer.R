# glimmer
# translate glmer style model formulas into equivalent formula alist suitable for map2stan
# just a stub for now

detectunivar <- function(name) {
  res <- grep(name, pattern = '^((?<![c(]).)*[[]',perl = TRUE)
  if(length(res) > 0)
    return(TRUE)
  else return(FALSE)
}
    
undot <- function( astring ) {
    astring <- gsub( "." , "_" , astring , fixed=TRUE )
    astring <- gsub( ":" , "_X_" , astring , fixed=TRUE )
    astring <- gsub( "(" , "" , astring , fixed=TRUE )
    astring <- gsub( ")" , "" , astring , fixed=TRUE )
    astring
}

concat <- function( ... ) {
  paste( ... , collapse="" , sep="" )
}

coerce_index <- function( x ) as.integer(as.factor(as.character(x)))

nobars <- function(term) {
  if (!('|' %in% all.names(term))) return(term)
  if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
  if (length(term) == 2) {
    nb <- nobars(term[[2]])
    if (is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars(term[[2]])
  nb3 <- nobars(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

findbars <- function(term) {
  if (is.name(term) || !is.language(term)) return(NULL)
  if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
  if (!is.call(term)) stop("term must be of class call")
  if (term[[1]] == as.name('|')) return(term)
  if (length(term) == 2) return(findbars(term[[2]]))
  c(findbars(term[[2]]), findbars(term[[3]]))
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
  if (is.name(term) || !is.language(term)) return(term)
  if (length(term) == 2) {
  term[[2]] <- subbars(term[[2]])
  return(term)
  }
  stopifnot(length(term) >= 3)
  if (is.call(term) && term[[1]] == as.name('|'))
  term[[1]] <- as.name('+')
  for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
  term
}

hasintercept <- function(term) {
    attr( terms(term) , "intercept" )==1
}

clean_name <- function(ranef_expr){
  return(undot(deparse(ranef_expr)))
}

expand_slash <- function(ranef_expr){
  ranef_name <- clean_name(ranef_expr[[3]])
  var <- list()
  if(regexec(pattern = '/', fixed = TRUE, text = ranef_name) == -1) {
    return(ranef_expr)
  } else {
    vars <- strsplit(ranef_name,split = '/',fixed=TRUE)[[1]]
    len_vars <- length(vars)
    if(len_vars > 2)
      stop('Only supports maximum of 1 "/"')
    interact_ranef_name <- paste(vars[1],vars[2],sep=':')
    ranef_expr[[3]] <- as.name(vars[1])
    interact_to_bars <- findbars(as.formula(paste('~ (1 | ',interact_ranef_name, ')', sep='')))
  }
  return(list(ranef_expr,interact_to_bars))
}

corr_re_test <- function(ranef_expr){
  if(length(ranef_expr[[2]]) == 1){
   return(TRUE) 
  } else if(deparse(ranef_expr[[2]][[2]]) != '0'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

build_ranefs <- function(var_list, corr_list, data){
  names_re <- sapply(var_list, function(i) clean_name(i[[3]]))
  names_re_list <- unique(names_re)
  num_re <- length(names_re_list)
  ranef <- sapply(1:num_re,function(i) setNames(vector(mode='list',length=2),c('corr','uncorr')),simplify=FALSE)
  ranef <- setNames(ranef,names_re_list)
  group_inds_list <- list()
  for(i in 1:length(var_list)){
    v <- var_list[[i]][[2]]
    name <- names_re[i]
    corr <- corr_list[[i]]
    v_class <- class(v)
    if (v_class =="numeric") {
      if (corr){
        ranef[[ name ]]$corr <- c("(Intercept)",ranef$corr[[name]])
      } else {
        ranef[[name]]$uncorr <- c("(Intercept)",ranef$uncorr[[name]])
      }
    } else {
        f <- as.formula( paste( "~-1+" , deparse( v ) , sep="" ) )
        if(corr){
          ranef[[ name ]]$corr <- c(ranef[[name]]$corr,colnames( model.matrix( f , data ) ))
        } else {
          ranef[[ name ]]$uncorr <- c(ranef[[name]]$uncorr,colnames( model.matrix( f , data ) ))
        }
    }
    group_inds_list[[name]] <- with(data, eval(var_list[[i]][[3]]))
  }
  return(list(ranef = ranef, group_inds_list = group_inds_list))
}

build_name <- function(re, prefix, group_var, num_group_vars) {
  aterm <- undot(re)
  newterm <- ""
  var_prefix <- prefix[2]
  if (num_group_vars>1) var_prefix <- concat(var_prefix, group_var, "_")
  if (aterm == "Intercept") {
      par_name <- concat(var_prefix, aterm)
      newterm <- concat(par_name, "[", group_var, "]")
  } else {
      par_name <- concat(var_prefix, aterm)
      newterm <- concat(par_name, "[", group_var, "]", "*", aterm)
  }
  return(list(newterm = newterm, par_name = par_name))
}

build_priors <- function(member, group_var, matt_trick) {
  prior_list <- list()
  gvar_name <- concat(member, "[", group_var, "]")
  prepend <- concat(member,'_')
  if(matt_trick){
    prior_list[[gvar_name]] <- concat("dnormstd(0,sigma_", prepend, group_var, ")")
  } else {
    prior_list[[gvar_name]] <- concat("dnorm(0,sigma_", prepend, group_var, ")")
  }
  prior_list[[concat("sigma_",prepend,group_var)]] <- concat("dcauchy(0,2)")
  return(prior_list)
}

apply_build_priors <- function(members, group_var, matt_trick) {
  prior_list <- NULL
  if(length(members) < 1)
    return(prior_list)
  for (i in 1:length(members))
    prior_list <- c(prior_list, build_priors(members[[i]], group_var, matt_trick))
  return(prior_list)
}

xparse_glimmer_formula <- function( formula , data ) {
    ## take a formula and parse into fixed effect and varying effect lists
    
    # find fixed effects list by deleting random effects and expanding
    f_nobars <- nobars(formula) ## Makes a single fixed effects formula
    # catch implied intercept error -- happens when right side of formula is only () blocks
    if (class(f_nobars)=="name" & length(f_nobars)==1) {
        f_nobars <- nobars(as.formula(paste(deparse(formula), "+ 1" )))
    }
    
    ## gets names of fixed  effects from colums of of data frame output by model.matrix
    
    des_frame_fixef <- as.data.frame(model.matrix( f_nobars , data ))
    fixef <- names(des_frame_fixef)  
    
    # convert to all fixed effects and build needed model matrix
    # used to be des_frame_fixef <- as.data.frame(model.matrix( subbars( formula ) , data ))
    outcome_name <- deparse(f_nobars[[2]])
    # des_frame_fixef <- cbind( data[[outcome_name]] , des_frame_fixef )
    # colnames(des_frame_fixef)[1] <- outcome_name
    outcome <- model.frame(f_nobars, data)[,1]
    if (class(outcome) == "matrix") {
        # fix outcome name
        outcome_name <- colnames(outcome)[1]
    }
    
    # check for any varying effects
    if (formula == nobars(formula)) {
        # no varying effects
        ranef <- list()
      var_mod <- list()
    } else {
        # find varying effects list
        var <- findbars(formula)
        var_mod <- sapply(var, function(i) expand_slash(i))
        var_mod <- unlist(var_mod, recursive=FALSE)
        corr_uncorr <- sapply(var_mod, function(i) corr_re_test(i))
        ranef <- build_ranefs(var_mod, corr_uncorr, data)
        inds <- do.call(cbind,ranef$group_inds_list)
        names(inds) <- names(ranef$group_inds_list)
        des_frame_fixef <- cbind(des_frame_fixef,inds)
    }
    
    # result sufficient information to build Stan model code
    list( y=outcome , yname=outcome_name , fixef=fixef , ranef=ranef$ranef , dat=des_frame_fixef, var = var_mod)
}


glimmer <- function( formula , data , family=gaussian , prefix=c("b_","v_") , default_prior="dnorm(0,10)" , matt_trick = TRUE, ... ) {
    
    undot <- function( astring ) {
        astring <- gsub( "." , "_" , astring , fixed=TRUE )
        astring <- gsub( ":" , "_X_" , astring , fixed=TRUE )
        astring <- gsub( "(" , "" , astring , fixed=TRUE )
        astring <- gsub( ")" , "" , astring , fixed=TRUE )
        astring
    }
    
    # convert family to text
    family.orig <- family
    if ( class(family)=="function" ) {
        family <- do.call(family,args=list())
    }
    link <- family$link
    family <- family$family
    
    # templates
    family_liks <- list(
        gaussian = "dnorm( mu , sigma )",
        binomial = "dbinom( size , p )",
        poisson = "dpois( lambda )"
    )
    lm_names <- list(
        gaussian = "mu",
        binomial = "p",
        poisson = "lambda"
    )
    link_names <- list(
        gaussian = "identity",
        binomial = "logit",
        poisson = "log"
    )
    
    # check input
    if ( class(formula)!="formula" ) stop( "Input must be a glmer-style formula." )
    if ( missing(data) ) stop( "Need data" )
    
    f <- formula
    flist <- alist()
    prior_list <- alist()
    
    # parse
    pf <- xparse_glimmer_formula( formula , data )
    pf$yname <- undot(pf$yname)
    
    # build likelihood
    # check for size variable in Binomial
    dtext <- family_liks[[family]]
    if ( family=="binomial" ) {
        if ( class(pf$y)=="matrix" ) {
            # cbind input
            pf$dat[[pf$yname]] <- as.integer(pf$y[,1])
            pf$dat[[concat(pf$yname,"_size")]] <- as.integer(apply(pf$y,1,sum))
            dtext <- concat("dbinom( ",concat(pf$yname,"_size")," , p )")
        } else {
            # bernoulli
            pf$dat[[pf$yname]] <- pf$y
            dtext <- concat("dbinom( 1 , p )")
        }
    } else {
        pf$dat[[pf$yname]] <- pf$y
    }
    flist[[1]] <- concat(as.character(pf$yname) ," ~ " , dtext)
    
    # build fixed linear model
    flm <- ""
    for ( i in 1:length(pf$fixef) ) {
        # for each term, add corresponding term to linear model text
        aterm <- undot(pf$fixef[i])
        newterm <- ""
        if ( aterm=="Intercept" ) {
            newterm <- aterm
            prior_list[[newterm]] <- default_prior
        } else {
            par_name <- concat( prefix[1] , aterm )
            newterm <- concat( par_name , "*" , aterm )
            prior_list[[par_name]] <- default_prior
        }
        if ( i > 1 ) flm <- concat( flm , " +\n        " )
        flm <- concat( flm , newterm )
    }
    
    vlm <- ""
    num_group_vars <- length(pf$ranef)
    if ( num_group_vars > 0 ) {
    for ( i in 1:num_group_vars ) {
      group_var <- undot(names(pf$ranef)[i])
      corr_members <- list()
      uncorr_members <- list()
      ranef <- c(pf$ranef[[i]]$corr,pf$ranef[[i]]$uncorr)
      num_corr <- length(pf$ranef[[i]]$corr)
      num_uncorr <- length(pf$ranef[[i]]$uncorr)
      for (j in 1:length(ranef)) {
        name_list <- build_name(ranef[[j]], prefix, group_var, num_group_vars)
        par_name <- name_list$par_name
        newterm <- name_list$newterm
        if (j <= num_corr){
          corr_members[[par_name]] <- par_name
        } else {
          uncorr_members[[par_name]] <- par_name
        }
        if ( i > 1 | j > 1 ) vlm <- concat( vlm , " +\n        " )
        vlm <- concat( vlm , newterm )
      }#j
        # add group prior
      if (length(corr_members)>1) {
        # multi_normal
        gvar_name <- concat( "c(" , paste(corr_members,collapse=",") , ")" , "[" , group_var , "]" )
        prior_list[[gvar_name]] <- concat( "dmvnormchol(0,sigma_" , group_var , ",L_Rho_" , group_var , ")" )
        prior_list[[concat("sigma_",group_var)]] <- concat( "dcauchy(0,2)" )
        prior_list[[concat("L_Rho_",group_var)]] <- concat( "dlkjcorrchol(2)" )
        prior_list_uncorr_members <- apply_build_priors(uncorr_members,group_var, matt_trick)
        prior_list <- c(prior_list, prior_list_uncorr_members)
      } else {
          prior_list_uncorr_members <- apply_build_priors(uncorr_members,group_var, matt_trick)
          if (length(corr_members) == 1)
            prior_list_corr_members <- build_priors(corr_members[[1]], group_var, matt_trick)
          prior_list <- c(prior_list, prior_list_corr_members, prior_list_uncorr_members)
      }
        # make sure grouping variables in dat
        # also ensure is an integer index
        pf$dat[[group_var]] <- coerce_index(pf$dat[,group_var])
        ## Used to be pf$dat[[group_var]] <- coerce_index( data[[group_var]] )
    }#i
    }# ranef processing
    
    # any special priors for likelihood function
    if ( family=="gaussian" ) {
        prior_list[["sigma"]] <- "dcauchy(0,2)"
    }
    
    # insert linear model
    lm_name <- lm_names[[family]]
    #link_func <- link_names[[family]]
    if ( vlm=="" )
        lm_txt <- concat( flm )
    else
        lm_txt <- concat( flm , " +\n        " , vlm )
    lm_left <- concat( link , "(" , lm_name , ")" )
    if ( link=="identity" ) lm_left <- lm_name
    flist[[2]] <- concat( lm_left , " <- " , lm_txt )
    
    # build priors
    for ( i in 1:length(prior_list) ) {
        pname <- names(prior_list)[i]
        p_txt <- prior_list[[i]]
        flist[[i+2]] <- concat( pname , " ~ " , p_txt )
    }
    
    # build formula text and parse into alist object
    flist_txt <- "alist(\n"
    for ( i in 1:length(flist) ) {
        flist_txt <- concat( flist_txt , "    " , flist[[i]] )
        if ( i < length(flist) ) flist_txt <- concat( flist_txt , ",\n" )
    }
    flist_txt <- concat( flist_txt , "\n)" )
    flist2 <- eval(parse(text=flist_txt))

    
    # clean variable names
    names(pf$dat) <- sapply( names(pf$dat) , undot )
    # remove Intercept from dat
    pf$dat[['Intercept']] <- NULL

    # result
    invisible(list(f=flist2,d=pf$dat))
}

stan_lmer <- function(fl , data , family=gaussian , prefix=c("b_","v_") , default_prior="dnorm(0,10)" , file_name=NULL, ... ){
  if(is.null(file_name)){
    warning("Generated stan code should go into a file with .stan suffix")
  }
  glim_result <- glimmer(fl , data , family=gaussian , prefix=c("b_","v_") , default_prior="dnorm(0,10)" , matt_trick = TRUE,... )
  stan_code <- map2stan(glim_result$f, glim_result$d, sample=FALSE)
  if(!is.null(file_name)){
    filewrite <- file(file_name)
    writeLines(stan_code$model, filewrite)
    close(filewrite)
  }
  invisible(list(stan_code = stan_code$model, data=stan_code$data))
}