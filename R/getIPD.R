#' Reconstruct individual patient data (IPD) from Kaplan-Meier (KM) survival curves
#'
#' This function returns reconstructed individual patient data (IPD) based on digitized data from Kaplan-Meier (KM) survival curves. This function takes as an input a vector of digitized censoring times and will incorporate them into the IPD exactly as specified.
#'
#' If the KM curve is horizontal at its end, either (1) the coordinates of the end of the curve can be included in t and S, OR (2) the x-coordinate of the end of the curve can be included in cens.t.
#'
#'
#' @param n The number at risk at time zero.
#' @param t A vector of the x coordinates of the bottoms of the drops in survival.
#' @param S A vector of the y coordinates corresponding to 't', between 0 and 1. These represent survival proportions at the times in t.
#' @param cens.t A vector of times at which patients were censored (optional).
#' @return A data.frame of reconstructed IPD with two columns; (1) t: survival time, (2) event: indicator of event, where 1=event and 0=censored.
#' @export
getIPD = function(n, t, S, cens.t = NA) {

  # Sort provided data
  t = sort(t);
  S = sort(S, decreasing=T);
  cens.t = sort(cens.t);

  # Correct negative times and survival proportions above 1
  t[t<0] = 0;
  S[S>1] = 1;

  # Add (t=0,S=1)
  t = c(0,t)
  S = c(1,S)

  # Create shell table for reconstructed IPD
  maxFU = max(t,cens.t, na.rm=T);
  IPD = data.frame(t = rep(maxFU, n), event = 0);

  # Keep track of the last patient for which IPD was reconstructed for
  p = 0;


  # Loop through entries in t
  for (i in 2:length(t)){

    if (is.na(cens.t[1])==F){
      # Determine the number censored between t[i-1] and t[i]
      n.cens = min(sum(cens.t>=t[i-1] & cens.t<t[i]), n-p);

      # Add censoring events to the IPD
      if (n.cens>0){
        IPD$t[(p+1):(p+n.cens)] = cens.t[cens.t>=t[i-1] & cens.t<t[i]][1:n.cens];
        IPD$event[(p+1):(p+n.cens)] = 0;
      }
      p = p+n.cens;
    }

    # Find which number of events at t[i] results in a survival proportion at t[i] closest to S[i]
    # Keep testing an increaing number of events until the survival proportion at t[i] falls below S[i]

    if (p<n){
    n.died = -1;
    diff=1;

    while (diff>0 & p+n.died+1<=n){

      diff_last = diff;
      n.died = n.died + 1;

      # Add events to the reconstructed IPD
      if (n.died>0){
        IPD$t[p+n.died] = t[i];
        IPD$event[p+n.died] = 1;
      }

      # Find resulting survival estimate at t[i]
      est.S = summary(survival::survfit(Surv(t, event=event) ~ 1, data = IPD),t=t[i])$surv;
      diff = est.S-S[i];
    }

    # Check whether n.died or n.died-1 resulted in a closer fit at t[i]
    # If n.died-1 resulted in a closer fit, remove one event from the reconstructed IPD
    if (n.died>0 & abs(diff_last)<abs(diff))  {
      IPD$t[p+n.died] = maxFU;
      IPD$event[p+n.died] = 0;
      n.died = n.died - 1;
    };

    # Update p
    p = p+n.died;
    }
  }

  if (is.na(cens.t[1])==F){
    # Add censored observations after last time in t
    n.cens = min(sum(cens.t>=max(t)),n-p);
    if (n.cens>0){
      IPD$t[(p+1):(p+n.cens)] = cens.t[cens.t>=max(t)][1:n.cens];
      IPD$event[(p+1):(p+n.cens)] = 0;
    }
  }

  return(IPD);
}
