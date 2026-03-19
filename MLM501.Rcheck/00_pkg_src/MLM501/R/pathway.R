# Evidence Pathway Engine: Test -> Treat Simulation
# Now with Health-Economic Impact Modeling

#' Simulate a clinical pathway from diagnosis to outcome.
#' @param prevalence Disease prevalence (0-1).
#' @param sensitivity Test sensitivity (0-1).
#' @param specificity Test specificity (0-1).
#' @param rr Relative risk reduction from treatment.
#' @param harm_rate Rate of treatment-related adverse events (default 0).
#' @param n Population size (default 1000).
#' @return A list with events_prevented, total_harms, and net_benefit.
#' @export
simulate_pathway <- function(prevalence, sensitivity, specificity, rr, harm_rate = 0, n = 1000) {
  tp <- n * prevalence * sensitivity
  fn <- n * prevalence * (1 - sensitivity)
  fp <- n * (1 - prevalence) * (1 - specificity)
  tn <- n * (1 - prevalence) * specificity
  treated <- tp + fp
  events_prevented <- tp * (1 - rr)
  total_harms <- treated * harm_rate
  ncb <- events_prevented - total_harms
  return(list(events_prevented = events_prevented, total_harms = total_harms, net_benefit = ncb))
}

#' Calculate Economic Impact and "Cost of Fragility"
#' @description Translates clinical outcomes into monetary values using QALY.
#' @param pathway_result The output from simulate_pathway()
#' @param cost_per_qaly Monetary value of one Quality-Adjusted Life Year (e.g., 30000 for NHS)
#' @param qaly_gain_per_event QALYs gained for one prevented event (e.g., 5 for preventing a stroke)
#' @param qaly_loss_per_harm QALYs lost per adverse event (e.g., 0.5)
#' @return A list containing the total economic value and the cost of fragility.
#' @export
calculate_economic_impact <- function(pathway_result, cost_per_qaly = 30000, qaly_gain_per_event = 5, qaly_loss_per_harm = 0.5) {
  
  # Value of events prevented
  total_qaly_gained <- pathway_result$events_prevented * qaly_gain_per_event
  value_of_benefit <- total_qaly_gained * cost_per_qaly
  
  # Cost of harms
  total_qaly_lost <- pathway_result$total_harms * qaly_loss_per_harm
  cost_of_harm <- total_qaly_lost * cost_per_qaly
  
  # Net Economic Value
  net_economic_value <- value_of_benefit - cost_of_harm
  
  return(list(
    net_economic_value = net_economic_value,
    value_of_benefit_formatted = paste0("GBP", format(round(value_of_benefit/1e6, 2), nsmall=2), "M"),
    cost_of_harm_formatted = paste0("GBP", format(round(cost_of_harm/1e3, 0)), "K")
  ))
}

# [plot_pathway_frontier remains the same]
