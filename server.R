
# This is the server logic for a Shiny web application.
#  can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)


require(shiny)

# This function is used in the solver function and has no independent usages
MacroModeleq <- function(t, y, parms)
{
  with(
    as.list(c(y,parms)), #lets us access variables and parameters stored in y and pars by name
    {
      
      #the ordinary differential equations
      dH = (a - b) * H - (alpha + delta) * P  #host population
      dW =  lambda * P - gamma * W - beta * W * H #free living parasite stage
      dP =  beta * W * H - (u + b + alpha) * P - alpha * (P^2)/H *((k + 1)/k) #adult parasite
      
      list(c(dH, dW, dP))
    }
  ) #close with statement
} #end function specifying the ODEs


#' Simulation of a compartmental macroparasite  model illustrating 
#'
#' @description  This model allows for the simulation of macroparasite - host dynamics
#' 
#'
#' @param H0 initial number of hosts
#' @param W0 initial number of free living stages
#' @param P0 initial number of adult parasites
#' @param tmax maximum simulation time, units of weeks
#' @param a instantaneous birth rate of host
#' @param b instantaneous death rate of host due to all causes except parasite
#' @param alpha instantaneous death rate of host due to parasite
#' @param delta instantaneous reduction in host fecundity
#' @param lambda instantaneous birth rate of parasite eggs
#' @param gamma instantaneous death rate of free-living parasite stages
#' @param beta instantaneous rate of ingestion of parasite infective stages
#' @param u instantenous death rate of adult parasites 
#' @param k parameter of the negative binomial distribution which measures inversely 
#' the degree of aggregation of parasites within host population
#' @return This function returns the simulation result as obtained from a call
#'   to the deSolve ode solver
#' @details A compartmental ID model with several states/compartments
#'   is simulated as a set of ordinary differential
#'   equations. The function returns the output from the odesolver as a matrix,
#'   with one column per compartment/variable. The first column is time.
#' @section Warning:
#'   This function does not perform any error checking. So if you try to do
#'   something nonsensical (e.g. any negative values or fractions > 1),
#'   the code will likely abort with an error message
#' @examples
#'   # To run the simulation with default parameters just call this function
#'   result <- simulate_MacroModeltransmission()
#'   # To choose parameter values other than the standard one, specify them e.g. like such
#'   result <- simulate_MacroModeltransmission(H0 = 100, W0 = 1e5,  tmax = 100)
#'   # You should then use the simulation result returned from the function, e.g. like this:
#'   plot(result[,1],result[,2],xlab='Time',ylab='Number Hosts',type='l')
#' @seealso The UI of the shiny app 'Macroparasite Model', which is part of this package, contains more details on the model
#' @author Ania Majewska
#' @references See e.g. the article "Regulation and stability of a free-living host-parasite system: Trichostrongylus tenuis in red grouse. II. Population models." by Dobson and Hudson 
#' for information on models of this type 
#' see the documentation for the deSolve package for details on ODE solvers
#' @export



simulate_MacroModeltransmission <- function(H0 = 50, W0 = 200, P0 = 200, tmax = 50, a = 1.1, b = 0.75,  alpha = 0.04, delta = 0.001, lambda = 11, gamma = 6 ,beta = 0.4, u = 1, k = 1)
{
  ############################################################
  Y0 = c(H = H0, W = W0, P = P0);  #combine initial conditions into a vector
  dt = min(0.1, tmax / 1000); #time step for which to get results back
  timevec = seq(0, tmax, dt); #vector of times for which solution is returned (not that internal timestep of the integrator is different)
  
  ############################################################
  #vector of parameters which is sent to the ODE function  
  pars=c(a = a, b = b,  alpha = alpha, delta = delta, lambda = lambda, gamma = gamma, beta = beta, u = u, k = k); 
  
  #this line runs the simulation, i.e. integrates the differential equations describing the infection process
  #the result is saved in the odeoutput matrix, with the 1st column the time, the remaining columns the values for the variables
  #returned in the order as specified in Y0 and the return from the solver function
  odeoutput = deSolve::lsoda(Y0, timevec, func = MacroModeleq, parms=pars, atol=1e-12, rtol=1e-12);
  
  return (odeoutput)
}

#' @title A helper function that takes result from the simulators and produces plots and text output
#'
#' @description This function generates plots and text to be displayed in the Shiny UI. 
#' This is a helper function. This function processes multiple simulation runs, supplied as a list
#' @param input the shiny app input structure
#' @param output the shiny app output structure
#' @param allres multiple runs of the simulation, supplied as a list. each entry in the list is expected to be a matrix 
#' @param varlist an optional list of vectors containing the names of variables 
#'    that should be plotted and processed independently. if not supplied, all variables will be shown in one plot.
#'    Also, final results for all variables will either be computed according to the groups of variable/compartment 
#'    names provided in varlist or if not provided, all variables will be used. 
#' @return output a list with plot, text and warn elements for display in a shiny UI
#' @details This function is called by the shiny server to produce output returned to the shiny UI
#' @author Andreas Handel
#' @export

generate_simoutput <- function(input,output,allres,varlist = NULL)
{
  
  #nplots contains the number of plots to be produced. 
  #Buy default it's a single plot
  #some apps ask for more than 1 plot, this is then sent in as a list to this function
  #check if user provided a list of variables to be processed separately
  nplots = 1; 
  if (!is.null(varlist))
  {
    nplots = length(varlist)
  }
  
  
  # Here, we use the result returned from the ode solver to produce the plot
  # the resulting plot is saved in the "plot" placeholder of the output variable
  output$plot <- renderPlot({
    input$submitBtn
    
    tmax = isolate(input$tmax)
    
    #make potentially more than one plot
    graphics::par(mfrow=c(nplots,1))
    
    #loop over variable sets for which to produce plots
    for (vn in 1:nplots)
    {
      #for multiple plots, names of variables to be plotted as passed in by varlist, otherwise naems are just all column names (minus time) 
      #for some reason only the <- operator works, not the = operator
      ifelse(nplots>1, varnames <- unlist(varlist[vn]), varnames <- colnames(allres()[[1]])[-1] )
      
      #if the app doesn't have an nreps setting, assign repetition = 1, otherwise use nreps setting
      nreps = ifelse(is.null(isolate(input$nreps)),1,isolate(input$nreps)) 
      
      #process first simulation to build plot
      res = allres()[[1]]      
      ymax = max(res[,varnames]) #each subplot gets a separate ymax
      tvec = res[,1]
      
      
      mycols=c("blue",'orange','red','green','black','magenta','cyan')
      #plot the 1st line
      graphics::plot(tvec,res[,varnames[1]],type="l",xlab="time",ylab="",col=mycols[1],lwd=1,log="",xlim=c(0,tmax),ylim=c(0,ymax),main="Time Series")
      
      if (length(varnames)>1) #plot additional lines if there is more than 1 variable to be plotted
      {
        for (nn in 2:length(varnames))
        {
          graphics::lines(tvec,res[,varnames[nn]],type="l",col=mycols[nn],lwd=1,lty=1)
        }
      }
      graphics::legend("right", varnames,col = mycols,lty=c(1),lwd=2)
      
      
      #loop over each additional simulation
      if (nreps>1)
      {
        #results are added to plot
        for (n1 in 2:nreps)
        {
          res = allres()[[n1]]      
          tvec = res[,1]
          graphics::lines(tvec,res[,varnames[1]],type="l",col=mycols[1],lwd=1,lty=1) #first variable for each new simulation
          if (length(varnames)>1) #plot additional lines if there is more than 1 variable to be plotted
          {
            for (nn in 2:length(varnames))
            {
              graphics::lines(tvec,res[,varnames[nn]],type="l",col=mycols[nn],lwd=1,lty=1)
            }
          } #done adding additional variables
        } #done additing lines from additional runs
      } #end loop over addition simulation replications
      
    } #end loop over individual plots
    
  } #finish render-plot statement
  , width = 'auto', height = 'auto'
  ) #end the output$plot function which produces the plot
  
  # Use the result "res" returned from the simulator to compute some text results
  # the text should be formatted as HTML and placed in the "text" placeholder of the UI
  output$text <- renderUI({
    
    
    #if the app doesn't have an nreps setting, assign repetition = 1, otherwise use nreps setting
    nreps = ifelse(is.null(isolate(input$nreps)),1,isolate(input$nreps)) 
    
    #process sets of variables independently
    alltext <- ""
    
    #if multiple plots are requested, text output for variables will be processed 
    #using the same variable groupings as for the plots
    for (vn in 1:nplots) 
    {    
      #for multiple plots, names of variables to be plotted as passed in by varlist, otherwise names are just all column names (minus time) 
      ifelse(nplots>1, varnames <- unlist(varlist[vn]), varnames <- colnames(allres()[[1]])[-1] )
      
      resfinal = rep(0,length(varnames)) 
      resmax = rep(0,length(varnames))
      resmin = rep(0,length(varnames))
      resfracfinal = rep(0,length(varnames)) 
      for (n1 in 1:nreps) #add all final values
      {
        currentsim = allres()[[n1]]
        nrows = nrow(currentsim) #number of entries in time-series matrix - can be different for every run
        currfinal = currentsim[nrows,varnames] #final number for each variable of interest
        #min and max for each variable
        if (length(varnames)>1)
        {
          resmax = resmax + apply(currentsim[,varnames],2,max);
          resmin = resmin + apply(currentsim[,varnames],2,min);
        }
        if (length(varnames)==1) #for a single variable, we have a vector and the apply function does not work
        {
          resmax = resmax + max(currentsim[,varnames]);
          resmin = resmin + min(currentsim[,varnames]);
        }
        
        resfinal = resfinal + currfinal #total numbers
        resfracfinal = resfracfinal + currfinal / sum(currfinal) #add up fractions
      }  
      resmax = resmax/nreps; #mean across simulations (for stochastic models)
      resmin = resmin/nreps; #mean across simulations (for stochastic models)
      resfinal = resfinal/nreps #mean for each variable
      resfracfinal = resfracfinal/nreps #mean for each variable
      
      
      for (nn in 1:length(varnames))
      {
        maxval = round(resmax[nn],2)
        minval = round(resmin[nn],2)
        numfinal = round(resfinal[nn], 2)
        fracfinal = round(resfracfinal[nn], 2)
        newtxt1 <- paste('Minimum and Maximum of ',varnames[nn],' during simulation: ',minval,' and ', maxval,sep='')
        newtxt2 <- paste('Number and Fraction of ',varnames[nn],' at end of simulation: ',numfinal,' and ',fracfinal,sep='')
        if (nn == 1) {txt <- paste(newtxt1, newtxt2, sep = "<br/>")}
        if (nn > 1) {txt <- paste(txt, newtxt1, newtxt2, sep = "<br/>")}
      }
      alltext <- paste(alltext, txt, sep = "<hr>" ) #add text blocks together
      
    } #finishes loop over sets of variables
    
    finaltxt <- '<hr> <i> For stochastic simulation scenarios, values shown are the mean over all simulations. </i>'
    resulttxt <- paste(alltext, finaltxt, sep = "")
    HTML(resulttxt)
  }) #end text output
  
  # At last, if we have any warnings or error from the simulator we can show them here
  # That text will be shown in red in the UI ("warn" placeholder will be used)
  output$warn <- renderUI({
    warntxt <- ""
    if(length(utils::data()$warns) == 0){
      
    }else{
      warntxt <- paste(warntxt, "Warnings:", sep = "<br/>")
      for (i in 1:length(utils::data()$warns)){
        warntxt <- paste(warntxt, utils::data()$warns[[i]], sep = "<br/>")
      }
    }
    HTML(warntxt)
  })
}

#the server-side function with the main functionality
#this function is wrapped inside the shiny server function below to allow return to main menu when window is closed
refresh <- function(input, output){
  
  # This reactive takes the input data and sends it over to the simulator
  # Then it will get the results back and return it as the "res" variable
  res <- reactive({
    input$submitBtn
    
    # Read all the input values from the UI
    H0 = isolate(input$H0);
    W0 = isolate(input$W0);
    P0 = isolate(input$P0);
    tmax = isolate(input$tmax);
    a = isolate(input$a);
    b = isolate(input$b);
    alpha = isolate(input$alpha);
    delta = isolate(input$delta);
    lambda = isolate(input$lambda);
    gamma = isolate(input$gamma);
    beta = isolate(input$beta);
    u = isolate(input$u);
    k = isolate(input$k);
    
    
    # Call the ODE solver with the given parameters
    result <- simulate_MacroModeltransmission(H = H0, W = W0, P = P0, tmax = tmax, a = a, b = b, alpha = alpha, delta = delta, lambda = lambda, gamma = gamma, beta = beta, u = u, k = k)
    
    return(list(result)) #this is returned as the res variable
  })
  
  #if we want certain variables plotted and reported separately, we can specify them manually as a list
  #if nothing is specified, all variables are plotted and reported at once
  varlist = list(c("H","W","P"))
  #function that takes result saved in res and produces output
  #output (plots, text, warnings) is stored in and modifies the global variable 'output'
  generate_simoutput(input,output,res,varlist=varlist)
} #ends the 'refresh' shiny server function that runs the simulation and returns output

#main shiny server function
server <- function(input, output, session) {
  
  # Waits for the Exit Button to be pressed to stop the app and return to main menu
  observeEvent(input$exitBtn, {
    input$exitBtn
    stopApp(returnValue = 0)
  })
  
  # This function is called to refresh the content of the Shiny App
  refresh(input, output)
  
  # Event handler to listen for the webpage and see when it closes.
  # Right after the window is closed, it will stop the app server and the main menu will
  # continue asking for inputs.
  session$onSessionEnded(function(){
    stopApp(returnValue = 0)
  })
} #ends the main shiny server function

