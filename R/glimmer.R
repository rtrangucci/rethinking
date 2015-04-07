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

expand_slash <- function(ranef_expr){

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
        var <- findbars( formula )
        ranef <- list()
				sapply(var,function(i) print(i))
				i <- 1
        while (i <= length(var)) {
          str_name <- deparse( var[[i]][[3]] )
					name <- undot(str_name)
					if(regexec(pattern = '/', fixed = TRUE, text = name) != -1) {
						vars <- strsplit(name,split = '/',fixed=TRUE)[[1]]
						len_vars <- length(vars)
						if(len_vars > 2)
							stop('Only supports maximum of 1 "/"')
						interact_name <- paste(vars[1],vars[2],sep=':')
						var[[i]][[3]] <- as.name(vars[1])
						interact_to_bars <- findbars(as.formula(paste('~ (1 | ',interact_name, ')', sep='')))
						var <- c(var,interact_to_bars)
						name <- vars[1]
					} else if(ifelse(length(var[[i]][[2]]) > 1, 
													 ifelse(deparse(var[[i]][[2]][[2]]) == '0', TRUE, FALSE), FALSE)){
					}
					des_frame_fixef[[name]] <- with(data, eval(var[[i]][[3]]))
					
					# parse formula
					v <- var[[i]][[2]]
					if ( class(v)=="numeric" ) {
							# just intercept
							ranef[[ name ]] <- c("(Intercept)",ranef[[name]])
					} else {
							# should be class "call" or "name"
							# need to convert formula to model matrix headers
							f <- as.formula( paste( "~-1+" , deparse( v ) , sep="" ) )
							ranef[[ name ]] <- c(ranef[[name]],colnames( model.matrix( f , data ) ))
					}
					i <- i + 1
        }
    }
    
    # result sufficient information to build Stan model code
    list( y=outcome , yname=outcome_name , fixef=fixef , ranef=ranef , dat=as.data.frame(des_frame_fixef), var = var)
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
    flist[[1]] <- concat( as.character(pf$yname) , " ~ " , dtext )
    
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
				print(group_var)
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
						print(par_name)
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