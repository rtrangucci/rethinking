# glimmer
# translate glmer style model formulas into equivalent formula alist suitable for map2stan
# just a stub for now

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
	num_re <- length(var_list)
	names_re <- sapply(var_list, function(i) clean_name(i[[3]]))
  ranef <- sapply(1:num_re,function(i) setNames(vector(mode='list',length=2),c('corr','uncorr')),simplify=FALSE)
	ranef <- setNames(ranef,names_re)
	group_inds_list <- list()
	for(i in 1:num_re){
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


glimmer <- function( formula , data , family=gaussian , prefix=c("b_","v_") , default_prior="dnorm(0,10)" , ... ) {
    
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
            pf$dat[[pf$yname]] <- pf$y[,1]
            pf$dat[[concat(pf$yname,"_size")]] <- apply(pf$y,1,sum)
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
        members <- list()
        for ( j in 1:length(pf$ranef[[i]]) ) {
            aterm <- undot(pf$ranef[[i]][j])
            newterm <- ""
            var_prefix <- prefix[2]
            if ( num_group_vars>1 ) var_prefix <- concat( var_prefix , group_var , "_" )
            if ( aterm=="Intercept" ) {
                par_name <- concat( var_prefix , aterm )
                newterm <- concat( par_name , "[" , group_var , "]" )
            } else {
                par_name <- concat( var_prefix , aterm )
                newterm <- concat( par_name , "[" , group_var , "]" , "*" , aterm )
            }
            members[[par_name]] <- par_name
            if ( i > 1 | j > 1 ) vlm <- concat( vlm , " +\n        " )
            vlm <- concat( vlm , newterm )
        }#j
        # add group prior
        if ( length(members)>1 ) {
            # multi_normal
            gvar_name <- concat( "c(" , paste(members,collapse=",") , ")" , "[" , group_var , "]" )
            prior_list[[gvar_name]] <- concat( "dmvnorm2(0,sigma_" , group_var , ",Rho_" , group_var , ")" )
            prior_list[[concat("sigma_",group_var)]] <- concat( "dcauchy(0,2)" )
            prior_list[[concat("Rho_",group_var)]] <- concat( "dlkjcorr(2)" )
        } else {
            # normal
            gvar_name <- concat( members[[1]] , "[" , group_var , "]" )
            prior_list[[gvar_name]] <- concat( "dnorm(0,sigma_" , group_var , ")" )
            prior_list[[concat("sigma_",group_var)]] <- concat( "dcauchy(0,2)" )
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
    cat(flist_txt)
    cat("\n")
    invisible(list(f=flist2,d=pf$dat))
    
}