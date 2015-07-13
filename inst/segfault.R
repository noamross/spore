library(spore)

parms = list(
  max_i = 100,
  is = 0:100,
  lambda = 0.001,
  lambda_ex = 0.2,
  alpha = 0.1,
  alpha_power = 1,
  mu = 0.01,
  r = 0.5,
  d = 0.01,
  K = 100,
  init_pop = 100,
  time_max = 100,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 1,
  micro_relax_steps = 0,
  project = FALSE,
  n_sims = 100,
  control_min = 0,
  control_max = 0,
  v = 50,
  c = 200,
  progress = TRUE,
  micro_record = 0
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


micro_state = c(100, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
shadow_state = c(97.57118, -174.8631)
time = 0
micro_state_c_stepto(micro_state = micro_state, parms = parms,control = 0, time = 0, timeto = 1000, record = 0, run = 1)
