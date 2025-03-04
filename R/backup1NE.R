# Setting up workspace ----
# If you are using Base R, this sets the file path to wherever this script is:
#setwd(getSrcDirectory()[1])
# On R-Studio, use:
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



ETBD_migrateSYM.NE = function(initialtree,
                           t = 100,
                           JmaxV = c(1000, 1000),
                           split = F,
                           bud = T,
                           siteN = 2,
                           DIST = "SRS",
                           psymp = c(0.1,0.1),
                           SADmarg = .1,
                           exparm = c(-0.7,-0.7),
                           NegExpEx = T,
                           ExpSpParm = .7,
                           ExpSpParm2 = -.1,
                           ExpSp = T,
                           splitparm = .5,
                           constantEX = .1,
                           migprob1 = 0,
                           migprob2 = 0,
                           exparm2 = c(.5,.5),
                           Asteroid = 40,
                           Asteroidimpact = c(-.2, -.2),
                           GROW = F,
                           initialsize = 100,
                           threshold = .1,
                           reggie = JmaxV
)


{{
  #monitor objects
  extincttotal = NULL
  trip2 <- NULL
  mig = list()
  extinctsp = list()
  symp = list()
  symptrip = list()
  exsp = list()
  trees = list()
  migrates = list()
  exty = list()
  richy <- c()
  con = c()
  yuppy = 1
  deadpool = JmaxV

  ##run these to run a time step individually for sim testing
  # ipa = 1
  # psymp = .15
  # pallo = .0
  # split = F
  # bud = T
  # allopatric = F
  # siteN = 2
  # probleave = .0
  # DIST = "SRS"
  # watchgrow = T
  # mig_percent = .3
  # SADmarg = .1
  # JmaxV = c(1000, 1000)
  # NegExpEx = T
  # exparm = -0.7
  # ExpSpParm = 2
  # ExpSp = F
  # splitparm  = .3
  # migprob = .4
  # constantEX = 0

  #### A few small function needed for the main function
  '%!in%' <- function(x,y)!('%in%'(x,y))

  myFun <- function(n = 5000000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(99999, n, TRUE)), sample(LETTERS, n, TRUE))
  }

  #yes the simulation will crash if you genrate more that 3 million species ...
  abcd <-myFun(50000)



  {

    ########### initialization for multiple sites ##########

    #initial species sizes
    if (siteN != 1){
      siteN = 1:siteN
      initialsize = initialsize
      names = as.factor(c(paste(
        "t", stringr::str_pad(siteN, 3, pad = "0"), sep = ""
      )))
      tree = ape::as.phylo(~ names)
      tree$edge.length = rep(1, length(siteN))

      ### starting tree
      ITtree <- tree

      ##starting species matrix
      matrix_list <- list()
      for (s in 1:length(siteN)) {
        q <- matrix(1)
        row.names(q) = tree$tip.label[s]
        matrix_list[[s]] = q
      }

      ##blank matrix to use later
      test1 <- list()
      for (s in 1:length(siteN)) {
        q <- matrix(
          c(1),
          nrow = 1,
          ncol = 1,
          byrow = FALSE,
          dimnames = list(c("start0"),
                          c("sizes"))
        )
        test1[[s]] = q
      }


      matrix_list0 <- matrix_list

      ###adding the sizes to the matrix
      for (o in 1:length(matrix_list0)) {
        matrix_list0[[o]] <- matrix_list0[[o]] + initialsize - 1
      }

      temptree <- ITtree
      matrix_list6 <- matrix_list0
      blank1 <- matrix(nrow = 0, ncol = 1)
      allo <- list()

      for (o in 1:length(siteN)) {
        allo[[o]] <- blank1
      }
    }
  }




  ########### initialization for ONE site ##########
  if (length(siteN) == 1){
    initialsize = initialsize
    site1 <- 1:2
    names1 = as.factor(c(paste(
      "t", stringr::str_pad(site1, 3, pad = "0"), sep = ""
    )))

    tree = ape::as.phylo(~ names1)
    tree$edge.length = rep(1, length(site1))
    ITtree <- tree

    matrix_list <- list()
    for (s in 1:length(siteN)) {
      q <- matrix(1)
      row.names(q) = tree$tip.label[s]
      matrix_list[[s]] = q
    }


    matrix_list0 <- matrix_list
    for (o in 1:length(matrix_list0)) {
      matrix_list0[[o]] <- matrix_list0[[o]] + initialsize - 1
    }


    test1 <- list()
    for (s in 1:length(siteN)) {
      q <- matrix(
        c(1),
        nrow = 1,
        ncol = 1,
        byrow = FALSE,
        dimnames = list(c("start0"),
                        c("sizes"))
      )
      test1[[s]] = q
    }


    temptree <- ITtree
    matrix_list6 <- matrix_list0
    blank1 <- matrix(nrow = 0, ncol = 1)
    allo <- list()

    for (o in 1:length(siteN)) {
      allo[[o]] <- blank1
    }
  }



  pine <- "(t001:1,t002:1);"



  print("version waterbuffalo")

  for (ipa in 1:t)

  {


    #start
    matrix_list0 <- matrix_list6


    matrix_list0
    ##deleting extinct species from matrix list 6 from previous step

    matrix_list05 <- DeleteExtinct(matrix_list0)
    matrix_list55 <- DeleteExtinct(matrix_list0)

    for( o in 1:length(matrix_list05)){
      matrix_list05[[o]] <- (na.exclude(matrix_list05[[o]]))
      attributes(matrix_list05[[o]])$na.action <- NULL
    }


    for( o in 1:length(matrix_list05)){
      if (length(matrix_list05[[o]]) == 0){
        print('site is extinct')
      }
    }





    ##### selecting species to migrate ######

    if (length(siteN) > 1) {
      dist <- makeLineDomain(length(siteN), migprob1)
      dist2 <- makeLineDomain(length(siteN), migprob2)
      dist[2,] <- dist2[2,]

      migratedata <- Migrate(matrix_list05, dist, 1, siteN)

      matrix_list1 <- migratedata$matrixlist

    } else {
      matrix_list1 <- matrix_list05
      migratedata <- c()

    }



    matrix_list05 <- DeleteExtinct(matrix_list1)

    for( o in 1:length(matrix_list05)){
      matrix_list05[[o]] <- (na.exclude(matrix_list05[[o]]))
      attributes(matrix_list05[[o]])$na.action <- NULL
    }

    matrix_list1 <- matrix_list05
    temptree <- pine


    ##adding species from allopatric speciaiton


    if (length(migratedata)>0) {

      Ne <- grow.newick(migratedata$allo, pine, abcd)
      Nallo.tree <- Ne$tree
      Atrip2 <- Ne$trip2
      abcd <- Ne$abcd


      ###adding allopatrically speciatig secies to matrix list

      matrix_list13 <- AlloSpec(matrix_list55,
                                migratedata$allo,
                                migratedata$old,
                                Atrip2,
                                splitparm,
                                siteN)
    } else {
      matrix_list13 <- matrix_list1
      Nallo.tree <- temptree
      Atrip2 <- c()
    }


    ##



    # Sympatric Speciation

    if (ExpSp) {

      stip = list()
      for (o in 1:length(matrix_list1)) {
        if (NA %!in% matrix_list1[[o]]) {
          # speciationp = ((matrix_list1[[o]][, 1])/JmaxV[o])^ExpSpParm
        #  speciationp = ((matrix_list1[[o]][, 1])/sum(unlist(matrix_list1[[o]])))^ExpSpParm
          speciationp = 1- exp(ExpSpParm2[o]*matrix_list1[[o]][, 1]^ExpSpParm[o])
          stip[[o]] <- speciationp
        }
      }



      #as logical...
      speciatinglog = list()
      for (o in 1:length(matrix_list1)) {
        if (!is.null(stip[[o]])) {
          if (NA %!in% (stip[[o]])) {
            splog = as.logical(rbinom(length(matrix_list1[[o]][,1]), 1, stip[[o]]))
            speciatinglog[[o]] <- splog
          } else {
            speciatinglog[[o]] <- matrix_list1[[o]]
          }
        }
      }
    } else {

      speciatinglog = list()
      for (o in 1:length(matrix_list1)) {
        spec = as.logical(rbinom(length(matrix_list1[[o]][,1]), 1, psymp[o]))   ##probability of sympatric speciation psymp
        ##probability of sympatric speciation psymp
        speciatinglog[[o]] = spec

      }
    }


    speciating = list()
    for (o in 1:length(siteN)) {
      i = 1
      if (length(matrix_list1[[o]]) != 0){
        spe = setdiff(rownames(matrix_list1[[o]])[speciatinglog[[o]]], extincttotal)
        speciating[[o]] = spe
      } else {
        speciating[[o]] = NA
      }
    }

    for (o in 1:length(siteN)) {
      if (NA %in% (speciating[[o]])) {
        speciating[[o]] <- character(0)
      }
    }

    symp_sp <- list()
    for (o in 1:length(speciating)) {
      speciatin <- subset(speciating[[o]], speciating[[o]] != "1")
      speciatin <- as.matrix(speciatin)
      symp_sp[[o]] <- speciatin
    }




    symp_sp <- DeleteDups(symp_sp)

    mag <- unlist(migratedata$allo)

    ##make migrated species unable to speciaiton sympatrically
    for (o in 1:length(symp_sp)) {
      for (l in 1:length(symp_sp[[o]])) {
        if (length(symp_sp[[o]][l]) > 0) {
          if (symp_sp[[o]][l] %in% mag) {
            symp_sp[[o]][l] <- NA
          }
        }
      }
    }

    for( o in 1:length(symp_sp)){
      symp_sp[[o]] <- (na.exclude(symp_sp[[o]]))
      attributes(symp_sp[[o]])$na.action <- NULL
    }

    matrix_list4 <- matrix_list13


    ###add the symp species onto the allotree
    Ne <- grow.newick(symp_sp, Nallo.tree, abcd)
    Nfull.tree <- Ne$tree
    trip2 <- Ne$trip2
    abcd <- Ne$abcd


    #temp hold all unique species
    temp <- unique(unlist(symp_sp))
    symptrip <- trip2

    #update the speciating tips
    trim <- trip2

    if (length(temp) > 0) {
      g <- 0
      tri <- list()
      triw <- list()
      for (o in 1:length(trim)) {
        tri[[g + 1]] <- trim[[o]][1]
        triw[[g + 1]] <- trim[[o]][2]
        g <- g + 1
      }


      splitparm <- runif(1, min = 0.1, max = 0.3)


      shell <- runif(1, min = 1, max = 100)

      if (shell >= 50) {
        bud = T
        split = F
      } else {
        bud = F
        split = T
      }


      ##for budding speciation one branch has same abundance and new branch has 10% of original
      if (bud) {
        i  <- 1
        fax <- list()
        for (o in 1:length(symp_sp)) {
          if (length(symp_sp[[o]]) != 0) {
            #update species names in sizes table
            m = match(symp_sp[[o]], row.names(matrix_list4[[o]]))
            for (k in 1:length(m)) {
              row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
              fax[i] <- (matrix_list4[[o]])[m][k]
              i <- i + 1
            }
          }
        }

        SV <- c()
        for (i in 1:length(matrix_list4)) {
          SV <- append(SV, length(matrix_list4[[i]]))
        }

        S <- SV

        J <- c()
        for ( i in 1:length(JmaxV)){
          J <- append(J, JmaxV[i] * (S[i] / (100 + S[i])))
        }

        oldJ <- c()
        for ( i in 1:length(matrix_list4)){
          old <- sum(matrix_list4[[i]])
          oldJ <- append(oldJ, old)
          deadpool <- J[[i]]-oldJ[[i]]
        }


        if (deadpool < 0){
          print("above")
          deadpool <- 30
        } else {
          print("below")
          deadpool <- deadpool}

        flop <- as.matrix(as.numeric(fax))
        random_percents <- runif(length(flop), min = 0, max = 1)  # Generates 10 random percentages

        farm <- c()
        for ( i in 1:length(flop)){
          clown <-  random_percents[i]*deadpool
          farm <- append(farm,clown)
          deadpool = deadpool - clown
        }

        flop <- as.matrix(as.numeric(farm))
      }

      ##for splitting speciation 10% is subtracted from original and new branch is 10% of original
      if (split) {
        i  <- 1
        fax <- list()
        for (o in 1:length(symp_sp)) {
          if (length(symp_sp[[o]]) != 0) {
            #update species names in sizes table
            m = match(symp_sp[[o]], row.names(matrix_list4[[o]]))
            for (k in 1:length(m)) {
              row.names(matrix_list4[[o]])[m][k] = paste(tri[i], sep = "")
              fax[i] <- (matrix_list4[[o]])[m][k]
              faax <- (matrix_list4[[o]])[m][k]
              faxtax <- as.numeric(faax) * splitparm
              matrix_list4[[o]][m][k] <-
                matrix_list4[[o]][m][k] - faxtax
              i <- i + 1
            }
          }
        }


      #10% of parent population abundance
      flop <- as.matrix(as.numeric(fax) * splitparm)
}
#




      ### pop is new species sizes and the new names
      pop <- symp_sp

      i <- 1
      for (o in 1:length(symp_sp)) {
        if (length(symp_sp[[o]]) >= 1) {
          for (k in 1:length(symp_sp[[o]])) {
            pop[[o]][k] <- flop[[i]]
            i <- i + 1
          }
        }
        pop[[o]] <- matrix(as.numeric(pop[[o]]))
      }

      i <- 1
      for (o in 1:length(symp_sp)) {
        if (length(symp_sp[[o]]) >= 1) {
          m = 1:length(symp_sp[[o]])
          for (k in 1:length(m)) {
            rownames(pop[[o]])[m][k] = paste(triw[i], sep = "")
            i <- i + 1
          }
        }
      }


      morto <- list()
      for (o in 1:length(siteN)) {
        for (h in 1:length(pop[[o]])) {
          mart <- rbind(matrix_list4[[o]], pop[[o]])
          morto[[o]] <- mart
        }
      }


      # Fix empty row names
      for (o in 1:length(siteN)) {
        for (k in length(morto[[o]]))
          if (length(rownames(morto[[o]])) < 1) {
            rownames(morto[[o]]) <- rownames(matrix_list4[[o]])
          }
      }
      matrix_list5 <- morto
    } else {
      matrix_list5 <- matrix_list4
    }


    if (NA %in% unlist(matrix_list5)) {
      message(
        "NA in matrixlist5: problem with adding sympatric species to matrix list",
        ipa
      )
    }

    ########updating survivors###########

    ##list of speciating species
    specspec <- unlist(symp_sp)
    specspec <- append(specspec, unlist(migratedata$allo))

    ##list of new species that have already grown in tree
    grownspec <- unlist(unique(symptrip))
    grownspec <- append(grownspec, unlist(unique(Atrip2)))
    allspec <-  unmatrixlist(matrix_list5)


    ##finding species that didn't speciate or grow or are extinct.
    sop = setdiff(allspec, specspec)
    sop2 = setdiff(sop, grownspec)
    sop3 = setdiff(sop2, extincttotal)
    ma <- test1
    ma[[1]] <- unique(sop3)


    #Updating survivors in tree
    Ntree <- Nfull.tree
    if (length(sop3 > 1)) {
      for (o in 1:length(siteN)) {
        for (k in 1:length(ma[[o]])) {
          if (ma[[o]][k] != 1) {
           # tree = surviveatx(tree, ma[[o]][k])
            Ntree = survive.NE(Ntree, ma[[o]][k])
          }
        }
      }
    }



preSAD <- matrix_list5



for( o in 1:length(matrix_list05)){
  matrix_list05[[o]] <- (na.exclude(matrix_list05[[o]]))
  attributes(matrix_list05[[o]])$na.action <- NULL
}


for( o in 1:length(matrix_list05)){
  if (length(matrix_list05[[o]]) == 0){
    print('site is extinct still')
  }
}


if (!GROW){

#print('still sad')

    ### RANKS ABUNDANCES AND DRAWS FROM SAD Fishers log series distribution
    if (DIST == "SRS") {
      if (length(unmatrixlist(matrix_list5)) > 5) {
        xx <- MakeSAD(matrix_list5, SADmarg, JmaxV)
        matrix_list5 <- xx
      }
    }


    if (DIST == "GEO") {
      if (length(unmatrixlist(matrix_list5)) > 5) {
        xx <- MakeGEO(matrix_list5, JmaxV)
        matrix_list5 <- xx

      }
    }

    if (DIST == "NORM") {
      if (length(unmatrixlist(matrix_list5)) > 5) {
        xx <- MakeLogNormalSAD(matrix_list5, SADmarg, JmaxV)
        matrix_list5 <- xx

      }
    }

    if (DIST == "NO") {
      if (length(unmatrixlist(matrix_list5)) > 5) {
        xx <- MakeNOSAD(matrix_list5, JmaxV)
        matrix_list5 <- xx

      }
    }




    if (NA %in% unlist(matrix_list5)) {
      message(
        "NA in matrixlist5: problem with the SAD rank setting",
        ipa
      )

    }


} else {

  xx <- NoForce(matrix_list5, JmaxV)
  matrix_list5 <- xx
}



library(TTR)


richy <- append(richy, length(unmatrixlist(matrix_list5)))

if (length(richy) > 50) {  # Ensure enough data points
  ema_values <- EMA(richy, n = 10)  # Smooth population data

  stable_count <- 0  # Track consecutive stable steps
  required_stable_steps <- 10  # How many steps need to be stable before stopping
  alpha <- threshold  # Set threshold as 1% of EMA

  for (t in (11:length(ema_values))) {  # Start from where EMA is valid
    if (!is.na(ema_values[t]) && !is.na(ema_values[t - 1])) {
      relative_change <- abs(ema_values[t] - ema_values[t - 1]) / ema_values[t]  # % change

      if (relative_change < alpha) {
        stable_count <- stable_count + 1
      } else {
        stable_count <- 0  # Reset if change is significant
      }

      if (stable_count >= required_stable_steps) {  # Ensure stability over time
        print(paste("Stopping at time step:", t))
        return(
          list(
            tree = Ntree,
            trees = trees,
            matrix_list = matrix_list6,
            mig = mig,
            migrates = migrates,
            exty = extinctsp,
            symp = symp,
            con = paste0("converged at ", ipa)
          )
        )
      }
    }
  }
}















for(o in 1:length(matrix_list5)) {
  if (length(matrix_list5[[o]]) == 1){
  if (is.na(matrix_list5[[o]])) {
    matrix_list5[[o]] <- preSAD[[o]]
  }
  }
}



if (ipa %in%  Asteroid) {
  JmaxV <- Asteroidimpact[yuppy]
  yuppy <- yuppy + 1

  print("asteroid hits")
} else {
  JmaxV <- reggie
}




    ####### extinction #########
    {

      # for (o in 1:length(matrix_list5)){
      #   for (k in 1:length(matrix_list5[[o]])){
      #     if (matrix_list5[[o]][k] == 0){
      #       matrix_list5[[o]][k] <- NA
      #     }
      #   }
      # }


    #  print(matrix_list5)

      #calculate extinction probability
      etip = list()
      for (o in 1:length(matrix_list5)) {
        if (NA %!in% matrix_list5[[o]]) {
          # extinctionp = 1 / (0.37 ^ 1.1) * exp(-1.1 * matrix_list5[[o]][, 1])
          if (NegExpEx) {
            #extinctionp = exp(exparm * matrix_list5[[o]][, 1])
            #extinctionp = exparm2*matrix_list5[[o]][, 1]^exparm
            extinctionp = 1- exp(exparm2[o]*matrix_list5[[o]][, 1]^exparm[o])
            etip[[o]] <- extinctionp
          } else {
            extinctionp = constantEX
            etip[[o]] <- extinctionp
          }
        }
      }


#print(extinctionp)


#
      # ##turning probabilities greater than 1 to 1
      # for (o in 1:length(matrix_list5)) {
      #   if (!is.null(etip[[o]])) {
      #     if (NA %!in% matrix_list5[[o]]) {
      #     for (i in 1:length(etip[[o]])) {
      #     #  if (!is.na(etip[[o]][i])) {
      #         if (etip[[o]][i] > 1) {
      #           etip[[o]][i] <- 1
      #         }
      #       }
      #     }
      #   }
      # }

      #as logical...
      etipp = list()
      for (o in 1:length(matrix_list5)) {
        if (!is.null(etip[[o]])) {
          if (NA %!in% (etip[[o]])) {
            extlog = as.logical(rbinom(length(matrix_list5[[o]][,1]), 1, etip[[o]]))
            etipp[[o]] <- extlog

          } else {
            etipp[[o]] <- matrix_list5[[o]]
          }
        }
      }


      #name and position of extinct
      extinct = list()
      for (o in 1:length(siteN)) {
        if (!is.null(etip[[o]][1])) {
          if (length(matrix_list5[[o]] > 0)) {
            ex <- row.names(matrix_list5[[o]])[etipp[[o]]]
            extinct[[o]] <- as.matrix(ex)
          }
        }
      }


      extinctx <- extinct

      ext <- c()
      for (o in 1:length(extinct)) {
        if (length(extinct[[o]]) > 0) {
          ext <- append(ext, (extinct[[o]]))
        }
      }

      matrix_list6 <- matrix_list5


      extincttotal <- c()

      for (o in 1:length(extinctx)) {
        if (length(extinctx) > 0) {
          if (length(extinctx[[o]]) > 0) {
            #prune the extinct from table
            m = match(extinctx[[o]], row.names(matrix_list5[[o]]))
            extincttotal = c(extincttotal,  matrix_list6[[o]][m, 1])
            matrix_list6[[o]][m, 1] = 0

          }
        }
      }


      deadpool = deadpool + sum(as.numeric(unlist(extincttotal)))

      # for (o in 1:length(matrix_list6)){
      #   if (NA %in% matrix_list6[[o]] ){
      #     matrix_list6[[o]] = c()
      #   }
      # }


      pine <- Ntree

      summat <- sum(na.omit(unlist(matrix_list6)))


      if (sum(summat) == 0){
        message('EVERYTHING IS DEAD')
        break
      }


      ## counting ten time steps
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
      B <- ipa/100
      if (is.wholenumber(B)){
        print(paste('relax everything is okay ...', ipa))
      }

    }


    #monitors of sizes and trees
    extinctsp[[ipa]] = extincttotal
    mig[[ipa]] = matrix_list6
    migrates[[ipa]] = migratedata$allo
     trees[[ipa]] = Ntree
     symp[[ipa]] = symp_sp

  }
  print(JmaxV)
  return(
    list(
      tree = Ntree,
      trees = trees,
      #final tree
      #all trees by timeslice
      matrix_list = matrix_list5,
      mig = mig,
      migrates = migrates,
      exty = extinctsp,
      symp = symp,
      con = c("not converged")
    )
  )

}
}


####### Extra testing #######
