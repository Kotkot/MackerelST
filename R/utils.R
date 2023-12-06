### Utility functions from lme4 packages 
  
  expandDoubleVerts <- function(term)
  {
    expandDoubleVert <- function(term) {
      frml <- formula(substitute(~x,list(x=term[[2]])))
      ## FIXME: do this without paste and deparse if possible!
      ## need term.labels not all.vars to capture interactions too:
      newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
      if(attr(terms(frml), "intercept")!=0)
        newtrms <- c("1", newtrms)
      
      as.formula(paste("~(",
                       paste(vapply(newtrms, function(trm)
                         paste0(trm, "|", deparse(term[[3]])), ""),
                         collapse=")+("), ")"))[[2]]
    }
    
    if (!is.name(term) && is.language(term)) {
      if (term[[1]] == as.name("(")) {
        term[[2]] <- expandDoubleVerts(term[[2]])
      }
      stopifnot(is.call(term))
      if (term[[1]] == as.name('||'))
        return( expandDoubleVert(term) )
      ## else :
      term[[2]] <- expandDoubleVerts(term[[2]])
      if (length(term) != 2) {
        if(length(term) == 3)
          term[[3]] <- expandDoubleVerts(term[[3]])
      }
    }
    term
  }
  
  findbars <- function(term)
  {
    ## Recursive function applied to individual terms
    fb <- function(term)
    {
      if (is.name(term) || !is.language(term)) return(NULL)
      if (term[[1]] == as.name("(")) return(fb(term[[2]]))
      stopifnot(is.call(term))
      if (term[[1]] == as.name('|')) return(term)
      if (length(term) == 2) return(fb(term[[2]]))
      c(fb(term[[2]]), fb(term[[3]]))
    }
    ## Expand any slashes in the grouping factors returned by fb
    expandSlash <- function(bb)
    {
      ## Create the interaction terms for nested effects
      makeInteraction <- function(x)
      {
        if (length(x) < 2) return(x)
        trm1 <- makeInteraction(x[[1]])
        trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
        list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
      }
      ## Return the list of '/'-separated terms
      slashTerms <- function(x)
      {
        if (!("/" %in% all.names(x))) return(x)
        if (x[[1]] != as.name("/"))
          stop("unparseable formula for grouping factor",call.=FALSE)
        list(slashTerms(x[[2]]), slashTerms(x[[3]]))
      }
      
      if (!is.list(bb))
        expandSlash(list(bb))
      else
        unlist(lapply(bb, function(x) {
          if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            ## lapply(unlist(...)) - unlist returns a flattened list
            lapply(unlist(makeInteraction(trms)),
                   function(trm) substitute(foo|bar, list(foo = x[[2]], bar = trm)))
          else x
        }))
    }## {expandSlash}
    
    modterm <- expandDoubleVerts(
      if(is(term, "formula")) term[[length(term)]] else term)
    expandSlash(fb(modterm))
  }
  
  nobars <- function(term) {
    nb <- nobars_(term)  ## call recursive version
    if (is(term,"formula") && length(term)==3 && is.symbol(nb)) {
      ## called with two-sided RE-only formula:
      ##    construct response~1 formula
      nb <- reformulate("1",response=deparse(nb))
    }
    ## called with one-sided RE-only formula, or RHS alone
    if (is.null(nb)) {
      nb <- if (is(term,"formula")) ~1 else 1
    }
    nb
  }
  
  nobars_ <- function(term)
  {
    if (!anyBars(term)) return(term)
    if (isBar(term)) return(NULL)
    if (isAnyArgBar(term)) return(NULL)
    if (length(term) == 2) {
      nb <- nobars_(term[[2]])
      if(is.null(nb)) return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobars_(term[[2]])
    nb3 <- nobars_(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }
  
  isBar <- function(term) {
    if(is.call(term)) {
      if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
        return(TRUE)
      }
    }
    FALSE
  }
  
  isAnyArgBar <- function(term) {
    if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
      for(i in seq_along(term)) {
        if(isBar(term[[i]])) return(TRUE)
      }
    }
    FALSE
  }
  
  anyBars <- function(term) {
    any(c('|','||') %in% all.names(term))
  }
  
  safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)
  abbrDeparse <- function(x, width=60) {
    r <- deparse(x, width)
    if(length(r) > 1) paste(r[1], "...") else r
  }
  
  barnames <- function(bars) vapply(bars, function(x) safeDeparse(x[[3]]), "")
  
  makeFac <- function(x,char.only=FALSE) {
    if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
  }
  
  reOnly <- function(f,response=FALSE) {
    response <- if (response && length(f)==3) f[[2]] else NULL
    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
                response=response)
  }
  
  RHSForm <- function(form,as.form=FALSE) {
    rhsf <- form[[length(form)]]
    if (as.form) reformulate(deparse(rhsf)) else rhsf
  }
  
  getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="", contrasts, sparse=FALSE) {
    ## fixed-effects model matrix X -
    ## remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))
    
    terms <- NULL ## make sure it's empty in case we don't set it
    
    nobs <- nrow(fr)
    
    ## check for empty fixed form
    ## need to ignore environments when checking!
    ##  ignore.environment= arg only works with closures
    idfun <- function(x,y) {
      environment(x) <- emptyenv()
      environment(y) <- emptyenv()
      return(identical(x,y))
    }
    
    if (idfun(RHSForm(fixedform, as.form=TRUE), ~ 0) ||
        idfun(RHSForm(fixedform, as.form=TRUE), ~ -1)) {
      X <- matrix(ncol=0, nrow=nobs)
      offset <- rep(0,nobs)
    } else {
      tt <- terms(fixedform)
      pv <- attr(mf$formula,"predvars")
      attr(tt, "predvars") <- fix_predvars(pv,tt)
      mf$formula <- tt
      terms_fixed <- terms(eval(mf,envir=environment(fixedform)))
      if (!sparse) {
        X <- model.matrix(drop.special2(fixedform), fr, contrasts)
      } else {
        X <- Matrix::sparse.model.matrix(drop.special2(fixedform), fr, contrasts)
        ## FIXME? ?sparse.model.matrix recommends MatrixModels::model.Matrix(*,sparse=TRUE)
        ##  (but we may not need it, and would add another dependency etc.)
      }
      ## will be 0-column matrix if fixed formula is empty
      offset <- rep(0,nobs)
      terms <- list(fixed=terms(terms_fixed))
      if (inForm(fixedform,quote(offset))) {
        ## hate to match offset terms with model frame names
        ##  via deparse, but since that what was presumably done
        ##  internally to get the model frame names in the first place ...
        for (o in extractForm(fixedform,quote(offset))) {
          offset_nm <- deparse(o)
          ## don't think this will happen, but ...
          if (length(offset_nm)>1) {
            stop("trouble reconstructing offset name")
          }
          offset <- offset + fr[[offset_nm]]
        }
      }
    }
    
    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    if (is.null(findbars(ranform))) {
      reTrms <- NULL
      Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
      ss <- integer(0)
    } else {
      
      ## FIXME: check whether predvars are carried along correctly in terms
      if (!ranOK) stop("no random effects allowed in ", type, " term")
      RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))
      
      mf$formula <- ranform
      reTrms <- mkReTrms(findbars(RHSForm(formula)), fr, reorder.terms=FALSE)
      
      ss <- splitForm(formula)
      ss <- unlist(ss$reTrmClasses)
      
      Z <- t(reTrms$Zt)   ## still sparse ...
    }
    
    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)
    
    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
    
    namedList(X, Z, reTrms, ss, terms, offset)
  }
  
  
  get_pars <- function(object, unlist = TRUE) {
    ee <- object$env
    x <- ee$last.par.best
    # work around built-in default to parList, which
    #  is bad if no random component
    if (length(ee$random)>0) x <- x[-ee$random]
    p <- ee$parList(x = x)
    # if (!unlist) return(p)
    # p <- unlist(p[names(p)!="b"])  ## drop primary RE
    # names(p) <- gsub("[0-9]+$","",names(p)) ## remove disambiguators
    p
  }
  
  ll_tweedie <- function(object, withheld_y, withheld_mu) {
    p <- stats::plogis(object$par[["thetaf"]]) + 1
    phi <- exp(object$par[["ln_phi"]])
    fishMod::dTweedie(y = withheld_y, mu = withheld_mu, p = p, phi = phi, LOG = TRUE)
  }
  
  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  
  
## Extract fixed effect variable names
  extract_variables <- function(formula){
    # extract everything in the parenthesis
    out <- unlist(strsplit(as.character(nobars(formula)), split='~')[3])
    # separate each variable and remove white space
    out <- unlist(strsplit(out, split='+', fixed=TRUE))
    out <- str_replace_all(out, fixed(" "), "")
    # remove intercept if existing
    if (TRUE %in% (c("0","-1","1") %in% out)) out <- out[-1]
    # extract everything in the parenthesis
    out <- unlist(sapply(out, function(x) ifelse(grepl("\\(", x), str_match(x, "(?<=\\().+?(?=\\))"), x)))
    # if the variable had a smoother, remove the "k" part
    out <- unlist(sapply(1:length(out), function(x) unlist(strsplit(out[x], split=',', fixed=TRUE))[1]))
    return(out)
  }
    

  ### Some functions 
  make_barrier_spde <- function(spde) {
    if ("spde_barrier" %in% names(spde)) {
      C0 <- spde$spde_barrier$C[[1]]
      C1 <- spde$spde_barrier$C[[2]]
      D0 <- spde$spde_barrier$D[[1]]
      D1 <- spde$spde_barrier$D[[2]]
      .I <- spde$spde_barrier$I
    } else {
      C0 <- rep(1, 2)
      C1 <- rep(1, 2)
      D0 <- Matrix::Matrix(0, 1, 1)
      D1 <- Matrix::Matrix(0, 1, 1)
      .I <- Matrix::Matrix(0, 1, 1)
    }
    list(C0 = C0, C1 = C1, D0 = D0, D1 = D1, I = .I)
  }  
  
  
###  some functions borrowed from the sdmTMB package
  make_year_i <- function(x) {
    x <- as.integer(as.factor(x))
    x - min(x)
  }    
  
  set_par_value <- function(opt, par) {
    as.numeric(opt[par == names(opt)])
  }   
  
  set_NA_map <- function(Map_phase2, var=c("SST", "CHL"), age=1) {
    if (length(var)==1) toNA <- grep(var, colnames(X))
    if (length(var)>1)  toNA <- sapply(var, function(x) grep(x, colnames(X)))
    Map_phase2$beta[toNA+ncol(X)*(age-1)] <- NA
    return(Map_phase2$beta)
  }    
  
  check_valid_factor_levels <- function(x, .name = "") {
    assert_that(is.factor(x),
                msg = sprintf("Random effect group column `%s` is not a factor.", .name))
    lev <- sort(levels(x))
    uni <- sort(unique(as.character(x)))
    assert_that(identical(lev, uni),
                msg = sprintf(
                  "Random effect group column `%s` has extra factor levels. Please remove them.", .name))
  }
  
  
## Utility function from VAST TMB
  Calc_Anisotropic_Mesh <-
    function(loc_x, loc_g, loc_i, Method, Extrapolation_List, anisotropic_mesh=NULL, fine_scale=FALSE, ...){
      
      #######################
      # Create the anisotropic SPDE mesh using 2D coordinates
      #######################
      
      # 2D coordinates SPDE
      if( fine_scale==FALSE ){
        if( is.null(anisotropic_mesh)){
          anisotropic_mesh = INLA::inla.mesh.create( loc_x, plot.delay=NULL, ...)
        }
      }else{
        loc_z = rbind( loc_x, loc_g, loc_i )
        outer_hull = INLA::inla.nonconvex.hull( as.matrix(loc_z), convex = -0.05, concave = -0.05)
        anisotropic_mesh = INLA::inla.mesh.create( loc_x, plot.delay=NULL, boundary=outer_hull, ...)
      }
      
      anisotropic_spde = INLA::inla.spde2.matern(anisotropic_mesh, alpha=2)
      
      # Exploring how to add projection-matrix A from knots to extrapolation-grid cells
      # Creating projection-matrix A is now done in `make_data`
      if( FALSE ){
        loc_g = as.matrix( Extrapolation_List$Data_Extrap[which(Extrapolation_List$Data_Extrap[,'Area_in_survey_km2']>0),c("E_km","N_km")] )
        outer_hull = INLA::inla.nonconvex.hull(loc_i, convex = -0.05, concave = -0.05)
        if( is.null(anisotropic_mesh)) anisotropic_mesh = INLA::inla.mesh.create( loc_x, plot.delay=NULL, boundary=outer_hull, ...)
        plot(anisotropic_mesh)
        A = INLA::inla.spde.make.A( anisotropic_mesh, loc_i )
        Check = apply( A, MARGIN=1, FUN=function(vec){sum(vec>0)})
        if( any(Check!=3) ) stop("Problem")
      }
      
      # Pre-processing in R for anisotropy
      Dset = 1:2
      # Triangle info
      TV = anisotropic_mesh$graph$tv       # Triangle to vertex indexing
      V0 = anisotropic_mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
      V1 = anisotropic_mesh$loc[TV[,2],Dset]
      V2 = anisotropic_mesh$loc[TV[,3],Dset]
      E0 = V2 - V1                      # E = edge for each triangle
      E1 = V0 - V2
      E2 = V1 - V0
      
      # Pre-processing for barriers
      # Barriers don't affect projection matrix A
      # Obtain polygon for water
      map_data = rnaturalearth::ne_countries( scale=switch("medium", "low"=110, "medium"=50, "high"=10, 50) )
      
      # Calculate centroid of each triangle in mesh and convert to SpatialPoints
      n_triangles = length(anisotropic_mesh$graph$tv[,1])
      posTri = matrix(NA, nrow=n_triangles, ncol=2)
      for(tri_index in 1:n_triangles){
        temp = anisotropic_mesh$loc[ anisotropic_mesh$graph$tv[tri_index,], ]
        posTri[tri_index,] = colMeans(temp)[c(1,2)]
      }
      posTri = sp::SpatialPoints(posTri, proj4string=sp::CRS(Extrapolation_List$projargs) )
      posTri = sp::spTransform(posTri, CRSobj=map_data@proj4string )
      
      # Calculate set of triangles barrier.triangles with centroid over land
      if( Method == "Barrier" ){
        anisotropic_mesh_triangles_over_land = unlist(sp::over(map_data, posTri, returnList=TRUE))
      }else{
        anisotropic_mesh_triangles_over_land = vector()
      }
      #
      #plot( x=posTri@coords[,1], y=posTri@coords[,2], col=ifelse(1:n_triangles%in%triangles_over_land,"black","red") )
      
      # Create Barrier object if requested
      # Don't do this unless necessary, because it sometimes throws an error
      #Diagnose issues:  assign("anisotropic_mesh", anisotropic_mesh, envir = .GlobalEnv)
      barrier_finite_elements = INLA:::inla.barrier.fem(mesh=anisotropic_mesh,
                                                        barrier.triangles=anisotropic_mesh_triangles_over_land)
      barrier_list = list(C0 = barrier_finite_elements$C[[1]],
                          C1 = barrier_finite_elements$C[[2]],
                          D0 = barrier_finite_elements$D[[1]],
                          D1 = barrier_finite_elements$D[[2]],
                          I = barrier_finite_elements$I )
      # sp::plot( INLA::inla.barrier.polygon(anisotropic_mesh, triangles_over_land) )
      
      # Calculate Areas 
      crossprod_fn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
      Tri_Area = rep(NA, nrow(E0))
      for(i in 1:length(Tri_Area)) Tri_Area[i] = crossprod_fn( E0[i,],E1[i,] )/2   # T = area of each triangle
      
      ################
      # Add the isotropic SPDE mesh for spherical or 2D projection, depending upon `Method` input
      ################
      
      # Mesh and SPDE for different inputs
      if(Method %in% c("Mesh","Grid","Stream_network","Barrier")){
        loc_isotropic_mesh = loc_x
        isotropic_mesh = anisotropic_mesh
      }
      if(Method %in% c("Spherical_mesh")){
        loc_isotropic_mesh = INLA::inla.mesh.map(loc_x, projection="longlat", inverse=TRUE) # Project from lat/long to mesh coordinates
        isotropic_mesh = INLA::inla.mesh.create( loc_isotropic_mesh, plot.delay=NULL, ...)
      }
      isotropic_spde = INLA::inla.spde2.matern(isotropic_mesh, alpha=2)
      
      ####################
      # Return stuff
      ####################
      #if( isotropic_mesh$n != anisotropic_mesh$n ) stop("Check `Calc_Anisotropic_Mesh` for problem")
      
      Return = list("loc_x"=loc_x, "loc_isotropic_mesh"=loc_isotropic_mesh, "isotropic_mesh"=isotropic_mesh,
                    "isotropic_spde"=isotropic_spde, "anisotropic_mesh"=anisotropic_mesh, "anisotropic_spde"=anisotropic_spde,
                    "Tri_Area"=Tri_Area, "TV"=TV, "E0"=E0, "E1"=E1, "E2"=E2,
                    "anisotropic_mesh_triangles_over_land"=anisotropic_mesh_triangles_over_land, "barrier_list"=barrier_list )
      return(Return)
    }
  

  
    