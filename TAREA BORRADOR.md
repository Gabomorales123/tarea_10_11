# =============================================================================
#  TAREA BORRADOR - Pipeline de análisis SedLine (R)            VERSIÓN 0.1
# -----------------------------------------------------------------------------
#  Tesis: Caracterización del monitoreo electroencefalográfico (SedLine)
#         durante la fase de transición de la ventilación mecánica en un
#         modelo porcino de SDRA (anidada en Fondecyt Regular 1261566).
#
#  !!!  ESTE SCRIPT TRABAJA SOBRE DATOS SINTÉTICOS  !!!
#  Los datos se generan a partir de parámetros publicados (Damiani et al.,
#  2022; Mirra et al., 2022a; Mirra et al., 2023b) y de SUPUESTOS marcados
#  como tales. Ninguna cifra, tabla o figura producida aquí es un resultado
#  real. Todos los archivos de salida llevan el prefijo "SINTETICO_".
#  Cuando existan datos reales, se reemplaza la SECCIÓN 2 (generación) por la
#  lectura de las exportaciones .csv del SedLine y de la CRF (SECCIÓN 3).
#
#  Correspondencia con los objetivos específicos (OE) de la propuesta:
#    OE1 Evolución horaria de PSI, SR, SEF95, EMG, ARTF (T0-T30) por grupo -> §5
#    OE2 Asociación intra-animal SedLine-esfuerzo (T18-T30)               -> §6-7
#    OE3 Comparación exploratoria entre grupos                            -> §7
#    OE4 Asociación SedLine - cambio de PdiTw (T30/T18)                   -> §8
#    OE5 Artefactos, contaminación por ECG y pérdida de datos             -> §4, §9
#
#  Alcance: estudio observacional y exploratorio (n = 21; 7 por grupo).
#  NO se ajustan modelos predictivos ni puntos de corte. Se reportan
#  estimaciones con IC 95 % y tamaños de efecto; los valores p son
#  descriptivos (sin ajuste confirmatorio por multiplicidad).
#  La TIE no forma parte de este análisis.
#
#  Reproducibilidad:
#    - Semilla fija (abajo).
#    - Entorno con renv (ejecutar UNA vez en la raíz del proyecto):
#        install.packages("renv"); renv::init(); renv::snapshot()
#    - Al final se guarda sessionInfo().
#
#  Estado: borrador NO ejecutado (R no está instalado en este equipo al
#  2026-11-10). Revisar después de la primera ejecución.
#  Sugerencia: renombrar este archivo a "pipeline_sedline.R".
# =============================================================================


# ---- 0. Configuración --------------------------------------------------------

SEMILLA <- 20261110
set.seed(SEMILLA)

paquetes <- c("dplyr", "tidyr", "ggplot2", "readr", "lme4", "lmerTest",
              "nlme", "rmcorr", "performance", "emmeans", "irr")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) {
  stop("Faltan paquetes: ", paste(faltan, collapse = ", "),
       "\nInstálalos con renv::install(c(\"", paste(faltan, collapse = "\", \""),
       "\")) y luego renv::snapshot().")
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(lme4); library(lmerTest); library(rmcorr)
  library(performance); library(emmeans); library(irr)
})
emm_options(lmer.df = "satterthwaite")

ETIQUETA <- "DATOS SINTÉTICOS - no corresponden a resultados reales"
PREFIJO  <- "SINTETICO_"

dir_out <- "salidas_sinteticas"
for (d in c("datos/sedline_crudo", "datos/derivados", "tablas", "figuras")) {
  dir.create(file.path(dir_out, d), recursive = TRUE, showWarnings = FALSE)
}
ruta <- function(...) file.path(dir_out, ...)

guardar_tabla <- function(df, nombre) {
  df <- dplyr::mutate(df, origen = "SINTETICO")
  readr::write_csv(df, ruta("tablas", paste0(PREFIJO, nombre, ".csv")))
  invisible(df)
}
guardar_figura <- function(p, nombre, ancho = 10, alto = 6) {
  ggsave(ruta("figuras", paste0(PREFIJO, nombre, ".png")), p,
         width = ancho, height = alto, dpi = 300)
}


# ---- 1. Parámetros -----------------------------------------------------------
# Fuente de cada parámetro indicada. "SUPUESTO" = valor plausible no tomado de
# un paper; debe revisarse antes de usar el pipeline como referencia.

P <- list(
  n_por_grupo     = 7,
  horas_total     = 30,          # T0-T30
  hora_transicion = 18,          # T18 = aleatorización
  dt_seg          = 2,           # SedLine guarda cada 2 s (propuesta, §3.3)

  # PTPmin (cmH2O·s·min^-1) - Damiani et al. (2022): nivel medio 185 ± 33,
  # nivel alto 601 ± 141. Se usan como metas de los grupos 2 y 3.
  ptp_protector = c(media = 185, de = 33),
  ptp_excesivo  = c(media = 601, de = 141),

  # PSI: meseta ~44,2 bajo propofol (Mirra et al., 2022a). La línea base en
  # sedación profunda se centra algo bajo la meseta (SUPUESTO).
  psi_meseta = 44.2,

  # Grado de "aligeramiento" de la sedación alcanzado en la transición
  # (0 = profunda, 1 = superficial). SUPUESTO.
  L_meta = c(pasivo = 0.05, protector = 0.55, excesivo = 0.90),
  tau_h  = 2.5,                  # constante de tiempo de la reducción (h). SUPUESTO

  # PdiTw T30/T18 - Damiani et al. (2022): esfuerzo alto -34 %, medio -7 %.
  # Grupo pasivo: SUPUESTO (-15 %).
  pditw_ratio = c(pasivo = 0.85, protector = 0.93, excesivo = 0.66),
  pditw_de    = 0.10,

  # Pérdida de datos del guardado automático. SUPUESTO.
  prob_hueco_por_hora = 0.08,
  dur_hueco_min       = c(2, 25),
  n_animales_ecg      = 3,       # animales con episodio de contaminación ECG

  # Control de calidad - UMBRALES PROVISIONALES: justificar con la guía
  # SedLine / Mirra et al. antes del análisis real.
  qc_artf_max     = 50,          # % ARTF sobre el cual la época se excluye
  qc_emg_alto     = 70,          # % EMG: se MARCA, no se excluye (ver §4)
  qc_ecg_sr_min   = 40,          # SR alto + PSI alto = posible contaminación ECG
  qc_ecg_psi_min  = 40,
  qc_gap_seg      = 10,          # salto de tiempo considerado hueco
  min_cobertura   = 0.50         # fracción mínima de épocas válidas por hora
)

GRUPOS <- c("pasivo", "protector", "excesivo")
epocas_por_hora <- 3600 / P$dt_seg


# ---- 2. Generación del dataset SINTÉTICO -------------------------------------
# (Reemplazar completa por la lectura de datos reales.)

n_total <- 3 * P$n_por_grupo
animales <- tibble(
  animal_id    = sprintf("C%02d", seq_len(n_total)),
  grupo        = factor(sample(rep(GRUPOS, each = P$n_por_grupo)), levels = GRUPOS),
  peso_kg      = round(rnorm(n_total, 30, 2.5), 1),
  psi_base     = rnorm(n_total, P$psi_meseta - 6, 4),  # efecto aleatorio animal
  sens         = rnorm(n_total, 1, 0.15),              # sensibilidad del PSI
  ecg_episodio = seq_len(n_total) %in% sample(seq_len(n_total), P$n_animales_ecg)
)

simular_animal <- function(a) {
  g    <- as.character(a$grupo)
  t_seg <- seq(0, P$horas_total * 3600 - P$dt_seg, by = P$dt_seg)
  t_h  <- t_seg / 3600
  n    <- length(t_seg)

  # Trayectoria latente de aligeramiento de la sedación
  L <- ifelse(t_h < P$hora_transicion, 0,
              P$L_meta[[g]] * (1 - exp(-(t_h - P$hora_transicion) / P$tau_h)))
  dev_h <- as.numeric(arima.sim(list(ar = 0.7), n = P$horas_total + 1, sd = 0.05))
  dev   <- approx(0:P$horas_total, dev_h, xout = t_h)$y
  L <- pmin(pmax(L + ifelse(t_h < P$hora_transicion, 0.3 * dev, dev), 0), 1)

  psi  <- a$psi_base + 30 * a$sens * L + rnorm(n, 0, 3)
  sr   <- pmax(0, rnorm(n, 6 * (1 - L) - 4, 4))
  sefl <- 14 + 1.5 * L + rnorm(n, 0, 2.5)        # SEF95 sin tendencia clara (Mirra 2022a)
  sefr <- sefl + rnorm(n, 0, 1)
  emg  <- 20 + 45 * L + rnorm(n, 0, 8)
  artf <- abs(rnorm(n, 0, 6))

  # Artefactos por manipulación: estimulación frénica en T0, T6 ... T30 (10 min)
  for (h in seq(0, P$horas_total, by = 6)) {
    idx <- t_h >= h & t_h < h + 10 / 60
    artf[idx] <- artf[idx] + runif(sum(idx), 40, 90)
  }

  # Contaminación por ECG en fase profunda (Mirra et al., 2023b): SR y PSI
  # falsamente elevados durante 30 min.
  ecg <- rep(FALSE, n)
  if (a$ecg_episodio) {
    ini <- runif(1, 2, 16)
    idx <- t_h >= ini & t_h < ini + 0.5
    ecg[idx] <- TRUE
    sr[idx]  <- runif(sum(idx), 45, 80)
    psi[idx] <- psi[idx] + runif(sum(idx), 15, 30)
  }

  completo <- tibble(
    animal_id = a$animal_id, t_seg = t_seg,
    PSI  = round(pmin(pmax(psi, 0), 100)),
    SR   = round(pmin(sr, 100)),
    SEFL = round(pmax(sefl, 0), 1),
    SEFR = round(pmax(sefr, 0), 1),
    EMG  = round(pmin(pmax(emg, 0), 100)),
    ARTF = round(pmin(artf, 100))
  )

  # Lectura manual en CRF a cada hora en punto (la pantalla muestra el valor
  # aunque falle el guardado). Error de transcripción y 3 % de faltantes.
  crf <- completo %>%
    filter(t_seg %% 3600 == 0) %>%
    transmute(animal_id, hora = t_seg / 3600,
              PSI_crf  = PSI + round(rnorm(n(), 0, 1.5)) +
                         ifelse(runif(n()) < 0.02, sample(c(-10, 10), n(), TRUE), 0),
              EMG_crf  = EMG + round(rnorm(n(), 0, 2)),
              SEFL_crf = round(SEFL + rnorm(n(), 0, 0.3), 1)) %>%
    mutate(across(ends_with("_crf"), ~ ifelse(runif(n()) < 0.03, NA, .x)))

  L_hora <- tibble(hora = floor(t_h), L = L) %>%
    group_by(hora) %>% summarise(L_media = mean(L), .groups = "drop")

  verdad_ecg <- tibble(animal_id = a$animal_id, t_seg = t_seg, ecg_verdadera = ecg)

  # Huecos del guardado automático
  exportado <- completo
  for (h in 0:(P$horas_total - 1)) {
    if (runif(1) < P$prob_hueco_por_hora) {
      ini <- h * 3600 + runif(1, 0, 3000)
      dur <- runif(1, P$dur_hueco_min[1], P$dur_hueco_min[2]) * 60
      exportado <- filter(exportado, !(t_seg >= ini & t_seg < ini + dur))
    }
  }

  list(sedline = exportado, crf = crf, L_hora = L_hora, verdad_ecg = verdad_ecg)
}

simular_esfuerzo <- function(a, L_hora) {
  g <- as.character(a$grupo)
  ptp_meta <- switch(g,
    pasivo    = 0,
    protector = max(50,  rnorm(1, P$ptp_protector[["media"]], P$ptp_protector[["de"]])),
    excesivo  = max(200, rnorm(1, P$ptp_excesivo[["media"]],  P$ptp_excesivo[["de"]])))
  n <- nrow(L_hora)
  en_trans <- L_hora$hora >= P$hora_transicion & g != "pasivo"

  ptp  <- ifelse(en_trans,
                 pmax(0, ptp_meta * L_hora$L_media / P$L_meta[[g]] + rnorm(n, 0, 0.1 * ptp_meta)),
                 0)
  dpes <- ifelse(ptp > 0, pmax(0, 0.035 * ptp + rnorm(n, 0, 1)), 0)
  dpdi <- ifelse(ptp > 0, pmax(0, 1.4 * dpes + rnorm(n, 0, 1)), 0)
  asin <- if (g == "excesivo") rnorm(n, 45, 12) else if (g == "protector") rnorm(n, 5, 3) else rep(0, n)
  asin <- ifelse(en_trans, pmin(pmax(asin, 0), 100), 0)

  tibble(animal_id = a$animal_id, hora = L_hora$hora,
         PTPmin = round(ptp, 1), dPes = round(dpes, 1), dPdi = round(dpdi, 1),
         asincronia_pct = round(asin, 1))
}

simular_pditw <- function(a) {
  g    <- as.character(a$grupo)
  base <- rnorm(1, 25, 4)                              # cmH2O, SUPUESTO
  horas <- c(0, 6, 12, 18, 24, 30)
  t18  <- base * (1 - 0.006 * 18)                      # leve caída en fase pasiva, SUPUESTO
  ratio <- P$pditw_ratio[[g]] * exp(rnorm(1, 0, P$pditw_de))
  val <- ifelse(horas <= 18, base * (1 - 0.006 * horas),
                ifelse(horas == 24, t18 * sqrt(ratio), t18 * ratio))
  tibble(animal_id = a$animal_id, hora = horas,
         PdiTw = round(val * exp(rnorm(length(horas), 0, 0.03)), 2))
}

sim <- lapply(seq_len(n_total), function(i) simular_animal(animales[i, ]))

# Exportación con formato "tipo dispositivo" (un archivo por animal)
for (i in seq_len(n_total)) {
  write_csv(sim[[i]]$sedline,
            ruta("datos/sedline_crudo", paste0(PREFIJO, animales$animal_id[i], ".csv")))
}
crf       <- bind_rows(lapply(sim, `[[`, "crf"))
verdad    <- bind_rows(lapply(sim, `[[`, "verdad_ecg"))
esfuerzo  <- bind_rows(lapply(seq_len(n_total), function(i)
               simular_esfuerzo(animales[i, ], sim[[i]]$L_hora)))
pditw     <- bind_rows(lapply(seq_len(n_total), function(i) simular_pditw(animales[i, ])))

write_csv(crf,      ruta("datos/derivados", paste0(PREFIJO, "crf_sedline_horaria.csv")))
write_csv(esfuerzo, ruta("datos/derivados", paste0(PREFIJO, "esfuerzo_horario.csv")))
write_csv(pditw,    ruta("datos/derivados", paste0(PREFIJO, "pditw.csv")))
write_csv(animales %>% select(animal_id, grupo, peso_kg),
          ruta("datos/derivados", paste0(PREFIJO, "animales.csv")))
rm(sim); invisible(gc())


# ---- 3. Lectura (punto de entrada para datos reales) -------------------------
# Con datos reales: apuntar a la carpeta de exportaciones y completar el mapeo
# de nombres de columnas del archivo de tendencias SedLine (verificar con la
# guía de referencia del equipo).

MAPEO_COLUMNAS <- c(PSI = "PSI", SR = "SR", SEFL = "SEFL", SEFR = "SEFR",
                    EMG = "EMG", ARTF = "ARTF")   # nombre_estandar = nombre_en_archivo

leer_sedline <- function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  d <- rename(d, any_of(MAPEO_COLUMNAS))
  d$animal_id <- sub(paste0("^", PREFIJO, "(C\\d+)\\.csv$"), "\\1", basename(f))
  d
}

archivos <- list.files(ruta("datos/sedline_crudo"), pattern = "\\.csv$", full.names = TRUE)
crudo <- bind_rows(lapply(archivos, leer_sedline)) %>%
  left_join(animales %>% select(animal_id, grupo), by = "animal_id")


# ---- 4. Control de calidad de la señal (OE5) ---------------------------------
# Reglas:
#  - ARTF > umbral       -> época excluida.
#  - Posible ECG         -> SR y PSI simultáneamente altos; excluida y
#                           marcada para verificación manual en el .edf.
#  - Fuera de rango / NA -> excluida.
#  - EMG alto            -> SOLO se marca. En los grupos con esfuerzo el EMG
#                           alto es parte del fenómeno estudiado; excluirlo
#                           sesgaría la comparación. Se evalúa en sensibilidad.

qc <- crudo %>%
  arrange(animal_id, t_seg) %>%
  group_by(animal_id) %>%
  mutate(hueco = c(0, diff(t_seg)) > P$qc_gap_seg) %>%
  ungroup() %>%
  mutate(
    flag_artf  = ARTF > P$qc_artf_max,
    flag_ecg   = SR >= P$qc_ecg_sr_min & PSI >= P$qc_ecg_psi_min,
    flag_rango = is.na(PSI) | PSI < 0 | PSI > 100,
    flag_emg   = EMG > P$qc_emg_alto,
    valido     = !(flag_artf | flag_ecg | flag_rango)
  )

esperadas <- P$horas_total * epocas_por_hora

qc_animal <- qc %>%
  group_by(animal_id, grupo) %>%
  summarise(
    n_registradas   = n(),
    pct_registrado  = 100 * n() / esperadas,
    n_huecos        = sum(hueco),
    pct_artf        = 100 * mean(flag_artf),
    pct_ecg         = 100 * mean(flag_ecg),
    pct_emg_alto    = 100 * mean(flag_emg),
    pct_analizable  = 100 * sum(valido) / esperadas,
    .groups = "drop")
guardar_tabla(qc_animal, "qc_por_animal")

qc_resumen <- qc_animal %>%
  summarise(across(starts_with("pct_"),
                   list(mediana = median, p25 = ~ quantile(.x, 0.25), p75 = ~ quantile(.x, 0.75))))
guardar_tabla(qc_resumen, "qc_resumen")

# Validación de la regla ECG (solo posible con datos sintéticos)
val_ecg <- qc %>%
  inner_join(verdad, by = c("animal_id", "t_seg")) %>%
  summarise(sensibilidad   = sum(flag_ecg & ecg_verdadera) / sum(ecg_verdadera),
            especificidad  = sum(!flag_ecg & !ecg_verdadera) / sum(!ecg_verdadera))
guardar_tabla(val_ecg, "validacion_regla_ecg")


# ---- 5. Agregación horaria y descriptivos por hora y grupo (OE1) -------------

VARS <- c("PSI", "SR", "SEFL", "SEFR", "EMG", "ARTF")

sedline_hora <- qc %>%
  mutate(hora = floor(t_seg / 3600)) %>%
  group_by(animal_id, grupo, hora) %>%
  summarise(
    n_registradas = n(),
    n_validas     = sum(valido),
    ARTF_todas    = mean(ARTF, na.rm = TRUE),
    across(all_of(VARS), ~ if (any(valido)) median(.x[valido], na.rm = TRUE) else NA_real_),
    .groups = "drop") %>%
  mutate(
    cobertura = n_validas / epocas_por_hora,
    across(all_of(VARS), ~ ifelse(cobertura >= P$min_cobertura, .x, NA_real_)),
    SEF95 = (SEFL + SEFR) / 2)
write_csv(sedline_hora, ruta("datos/derivados", paste0(PREFIJO, "sedline_horaria.csv")))

desc <- sedline_hora %>%
  pivot_longer(c(PSI, SR, SEF95, EMG, ARTF_todas), names_to = "variable", values_to = "valor") %>%
  group_by(grupo, hora, variable) %>%
  summarise(n       = sum(!is.na(valor)),
            mediana = median(valor, na.rm = TRUE),
            p25     = quantile(valor, 0.25, na.rm = TRUE),
            p75     = quantile(valor, 0.75, na.rm = TRUE),
            .groups = "drop")
guardar_tabla(desc, "descriptivos_hora_grupo")

p_desc <- ggplot(desc, aes(hora, mediana, colour = grupo, fill = grupo)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = P$hora_transicion, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Evolución horaria de las variables SedLine (mediana y RIC)",
       subtitle = ETIQUETA, x = "Hora (T0-T30)", y = NULL,
       caption = "Línea discontinua: T18, inicio de la fase de transición") +
  theme_minimal()
guardar_figura(p_desc, "fig1_evolucion_horaria")


# ---- 6. Correlación de medidas repetidas (OE2) -------------------------------
# rmcorr requiere variación intra-animal: el grupo pasivo tiene PTPmin = 0
# constante, por lo que se analiza en los grupos protector y excesivo.

trans <- sedline_hora %>%
  filter(hora >= P$hora_transicion, hora < P$horas_total) %>%
  inner_join(esfuerzo, by = c("animal_id", "hora")) %>%
  mutate(hora_c = hora - P$hora_transicion)

d_rm <- trans %>% filter(grupo != "pasivo") %>%
  mutate(animal_id = factor(animal_id)) %>% as.data.frame()

pares <- tribble(
  ~x,               ~y,
  "PTPmin",         "PSI",
  "dPdi",           "PSI",
  "dPes",           "PSI",
  "asincronia_pct", "PSI",
  "PTPmin",         "SEF95",
  "PTPmin",         "EMG",
  "PTPmin",         "SR"
)

res_rm <- bind_rows(lapply(seq_len(nrow(pares)), function(i) {
  d <- d_rm[complete.cases(d_rm[, c(pares$x[i], pares$y[i])]), ]
  r <- tryCatch(rmcorr(participant = "animal_id", measure1 = pares$x[i],
                       measure2 = pares$y[i], dataset = d),
                error = function(e) NULL)   # p. ej., SR sin variación intra-animal
  if (is.null(r)) return(tibble(x = pares$x[i], y = pares$y[i], r_rm = NA_real_,
                                n_animales = n_distinct(d$animal_id), n_obs = nrow(d)))
  tibble(x = pares$x[i], y = pares$y[i], r_rm = r$r,
         IC95_inf = r$CI[1], IC95_sup = r$CI[2], gl = r$df, p = r$p,
         n_animales = n_distinct(d$animal_id), n_obs = nrow(d))
}))
guardar_tabla(res_rm, "rmcorr_sedline_esfuerzo")

p_rm <- ggplot(d_rm, aes(PTPmin, PSI, colour = animal_id)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  facet_wrap(~ grupo, scales = "free_x") +
  labs(title = "PSI vs PTPmin dentro de cada animal (T18-T30)",
       subtitle = ETIQUETA, x = expression(PTP[min]~(cmH[2]*O%.%s%.%min^-1)), y = "PSI") +
  theme_minimal() + theme(legend.position = "none")
guardar_figura(p_rm, "fig2_rmcorr_psi_ptp")


# ---- 7. Modelos lineales mixtos (OE2-OE3) ------------------------------------
# Modelo principal: variable ~ grupo * hora + (1 | animal), T18-T30.
# Referencia = grupo protector. Sensibilidad: autocorrelación AR(1) (nlme).

trans_m <- trans %>% mutate(grupo = relevel(factor(grupo, levels = GRUPOS), ref = "protector"))

resumir_lmm <- function(m, nombre) {
  co <- summary(m)$coefficients
  ci <- confint(m, parm = "beta_", method = "Wald")
  r2 <- tryCatch(performance::r2_nakagawa(m), error = function(e) NULL)
  ic <- tryCatch(performance::icc(m),         error = function(e) NULL)
  tibble(modelo = nombre, termino = rownames(co),
         estimacion = co[, "Estimate"], EE = co[, "Std. Error"],
         IC95_inf = ci[rownames(co), 1], IC95_sup = ci[rownames(co), 2],
         p = co[, "Pr(>|t|)"],
         R2_marginal    = if (is.null(r2)) NA_real_ else r2$R2_marginal,
         R2_condicional = if (is.null(r2)) NA_real_ else r2$R2_conditional,
         ICC            = if (is.null(ic)) NA_real_ else ic$ICC_adjusted,
         singular       = isSingular(m))
}

contrastes_grupo <- function(m, nombre, horas = c(2, 11)) {
  sd_total <- sqrt(sum(as.data.frame(VarCorr(m))$vcov))
  emm <- emmeans(m, ~ grupo | hora_c, at = list(hora_c = horas))
  ctr <- summary(contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "none"),
                 infer = TRUE)
  as_tibble(ctr) %>%
    mutate(modelo = nombre,
           d_estandarizada = estimate / sd_total,
           d_IC95_inf = lower.CL / sd_total,
           d_IC95_sup = upper.CL / sd_total)
}

vars_lmm <- c("PSI", "SEF95", "EMG")
modelos  <- list()
tab_lmm  <- list()
tab_ctr  <- list()
for (v in vars_lmm) {
  f <- as.formula(paste(v, "~ grupo * hora_c + (1 | animal_id)"))
  m <- lmer(f, data = trans_m, REML = TRUE)
  modelos[[v]] <- m
  tab_lmm[[v]] <- resumir_lmm(m, v)
  tab_ctr[[v]] <- contrastes_grupo(m, v)
}
guardar_tabla(bind_rows(tab_lmm), "lmm_efectos_fijos")
guardar_tabla(bind_rows(tab_ctr), "lmm_contrastes_vs_protector")

# Sensibilidad AR(1) para PSI
d_ar <- trans_m %>% filter(!is.na(PSI))
m_lme  <- nlme::lme(PSI ~ grupo * hora_c, random = ~ 1 | animal_id,
                    data = d_ar, method = "REML")
m_ar1  <- nlme::lme(PSI ~ grupo * hora_c, random = ~ 1 | animal_id,
                    correlation = nlme::corAR1(form = ~ hora_c | animal_id),
                    data = d_ar, method = "REML")
sens_ar1 <- tibble(modelo = c("sin AR1", "con AR1"),
                   AIC = c(AIC(m_lme), AIC(m_ar1)),
                   phi = c(NA, coef(m_ar1$modelStruct$corStruct, unconstrained = FALSE)))
guardar_tabla(sens_ar1, "sensibilidad_ar1_psi")

# Descomposición intra/inter-animal del esfuerzo (grupos con esfuerzo)
wb <- trans_m %>% filter(grupo != "pasivo") %>%
  group_by(animal_id) %>%
  mutate(PTP_inter = mean(PTPmin, na.rm = TRUE) / 100,
         PTP_intra = (PTPmin - mean(PTPmin, na.rm = TRUE)) / 100) %>%
  ungroup()
m_wb <- lmer(PSI ~ PTP_intra + PTP_inter + hora_c + (1 | animal_id), data = wb)
guardar_tabla(resumir_lmm(m_wb, "PSI ~ PTPmin intra/inter (por 100 unidades)"),
              "lmm_psi_ptp_intra_inter")


# ---- 8. PdiTw (T30/T18) y SedLine (OE4) - exploratorio, n = 21 ---------------

pdi <- pditw %>%
  filter(hora %in% c(18, 30)) %>%
  pivot_wider(names_from = hora, values_from = PdiTw, names_prefix = "T") %>%
  mutate(ratio = T30 / T18, log_ratio = log(ratio))

expo <- trans %>%
  group_by(animal_id) %>%
  summarise(PSI_medio = mean(PSI, na.rm = TRUE),
            EMG_medio = mean(EMG, na.rm = TRUE),
            PTP_medio = mean(PTPmin, na.rm = TRUE), .groups = "drop")

d_pdi <- pdi %>% inner_join(expo, by = "animal_id") %>%
  left_join(animales %>% select(animal_id, grupo), by = "animal_id")
write_csv(d_pdi, ruta("datos/derivados", paste0(PREFIJO, "pditw_exposicion.csv")))

spearman_boot <- function(x, y, B = 2000) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  rho <- cor(x, y, method = "spearman")
  bs <- replicate(B, {
    i <- sample(seq_along(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  q <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  tibble(rho = rho, IC95_inf = q[[1]], IC95_sup = q[[2]], n = length(x))
}

res_pdi <- bind_rows(
  spearman_boot(d_pdi$PSI_medio, d_pdi$ratio) %>% mutate(x = "PSI_medio"),
  spearman_boot(d_pdi$EMG_medio, d_pdi$ratio) %>% mutate(x = "EMG_medio")
) %>% mutate(y = "PdiTw T30/T18")
guardar_tabla(res_pdi, "pditw_spearman")

m_pdi <- lm(log_ratio ~ PSI_medio, data = d_pdi)
guardar_tabla(tibble(termino = names(coef(m_pdi)), estimacion = coef(m_pdi),
                     IC95_inf = confint(m_pdi)[, 1], IC95_sup = confint(m_pdi)[, 2],
                     R2 = summary(m_pdi)$r.squared),
              "pditw_lm_logratio")

guardar_tabla(d_pdi %>% group_by(grupo) %>%
                summarise(n = n(), mediana_ratio = median(ratio),
                          p25 = quantile(ratio, 0.25), p75 = quantile(ratio, 0.75)),
              "pditw_por_grupo")

p_pdi <- ggplot(d_pdi, aes(PSI_medio, ratio, colour = grupo)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(size = 2.5) +
  labs(title = "Cambio de PdiTw (T30/T18) según PSI medio en la transición",
       subtitle = ETIQUETA, x = "PSI medio T18-T30", y = "PdiTw T30 / T18") +
  theme_minimal()
guardar_figura(p_pdi, "fig3_pditw_psi", ancho = 7, alto = 5)


# ---- 9. Concordancia registro automático vs CRF manual (OE5) -----------------
# Valor automático = mediana de épocas válidas en +/- 60 s de la hora en punto.

auto_hora <- qc %>%
  filter(valido) %>%
  mutate(hora = round(t_seg / 3600)) %>%
  filter(abs(t_seg - hora * 3600) <= 60) %>%
  group_by(animal_id, hora) %>%
  summarise(PSI_auto = median(PSI), EMG_auto = median(EMG),
            SEFL_auto = median(SEFL), .groups = "drop")

conc <- crf %>% left_join(auto_hora, by = c("animal_id", "hora"))

concordancia <- function(v) {
  d <- conc %>% select(auto = paste0(v, "_auto"), manual = paste0(v, "_crf")) %>%
    filter(!is.na(auto), !is.na(manual))
  ic <- irr::icc(as.matrix(d), model = "twoway", type = "agreement", unit = "single")
  dif <- d$manual - d$auto
  tibble(variable = v, n_pares = nrow(d),
         ICC_A1 = ic$value, ICC_IC95_inf = ic$lbound, ICC_IC95_sup = ic$ubound,
         sesgo = mean(dif), LdA_inf = mean(dif) - 1.96 * sd(dif),
         LdA_sup = mean(dif) + 1.96 * sd(dif))
}
# Nota: los límites de acuerdo ignoran las mediciones repetidas por animal;
# con datos reales usar la corrección de Bland y Altman para medidas repetidas.
guardar_tabla(bind_rows(lapply(c("PSI", "EMG", "SEFL"), concordancia)),
              "concordancia_auto_crf")

rescate <- conc %>%
  summarise(horas_con_crf        = sum(!is.na(PSI_crf)),
            horas_sin_auto       = sum(is.na(PSI_auto)),
            horas_rescatadas_crf = sum(!is.na(PSI_crf) & is.na(PSI_auto)))
guardar_tabla(rescate, "rescate_datos_por_crf")

d_ba <- conc %>% filter(!is.na(PSI_auto), !is.na(PSI_crf)) %>%
  mutate(media = (PSI_auto + PSI_crf) / 2, dif = PSI_crf - PSI_auto)
lda <- mean(d_ba$dif) + c(-1.96, 0, 1.96) * c(sd(d_ba$dif), 0, sd(d_ba$dif))
p_ba <- ggplot(d_ba, aes(media, dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = lda, linetype = c("dashed", "solid", "dashed")) +
  labs(title = "Bland-Altman: PSI manual (CRF) vs automático",
       subtitle = ETIQUETA, x = "Media de ambos registros", y = "CRF - automático") +
  theme_minimal()
guardar_figura(p_ba, "fig4_bland_altman_psi", ancho = 7, alto = 5)


# ---- 10. Registro del entorno -------------------------------------------------

writeLines(c(paste("Semilla:", SEMILLA), paste("Fecha:", Sys.time()), ETIQUETA, "",
             capture.output(sessionInfo())),
           ruta(paste0(PREFIJO, "sessionInfo.txt")))
message("Pipeline terminado. Salidas en: ", normalizePath(dir_out), "\n", ETIQUETA)
