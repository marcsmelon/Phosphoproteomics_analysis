################
# VOLCANO PLOT KSEA
##############

library(ggplot2)
library(dplyr)
library(ggrepel)

# ==============================================================================
# 1. CONFIGURACIÓN Y DATOS
# ==============================================================================
COL_DCM     <- "#F5CD6A"  # Oro (UP)
COL_HEALTHY <- "#3A2044"  # Morado (DOWN)
COL_NS      <- "#EAEAEA"  # Gris (No significativo)

# Dimensiones exactas (169.05 x 174.63 pt)
w_in <- 169.0586 / 72
h_in <- 174.6366 / 72

df <- read.csv("ksea_results.csv") %>%
  mutate(
    neg_log10_p = -log10(p_value),
    group = case_when(
      p_value < 0.05 & z_score > 0 ~ "UP",
      p_value < 0.05 & z_score < 0 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

# Conteos
n_up   <- sum(df$group == "UP")
n_down <- sum(df$group == "DOWN")

# Etiquetas: solo las 3 más significativas por grupo
top_labels <- df %>%
  filter(group != "NS") %>%
  group_by(group) %>%
  slice_min(order_by = p_value, n = 3) %>%
  ungroup()

# ==============================================================================
# 2. GRÁFICO (PUNTOS GRANDES, SIN BORDES)
# ==============================================================================
p <- ggplot(df, aes(x = z_score, y = neg_log10_p)) +
  # size = 1.5 hace los puntos más grandes y visibles
  # stroke = 0 asegura que no haya contorno
  geom_point(aes(color = group), alpha = 0.9, size = 1.5, stroke = 0) + 
  scale_color_manual(values = c("UP" = COL_DCM, "DOWN" = COL_HEALTHY, "NS" = COL_NS)) +
  
  # Líneas de referencia sutiles
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.2) +
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed", linewidth = 0.2) +
  
  # Anotaciones de conteo (n=)
  annotate("text", x = min(df$z_score)*0.8, y = max(df$neg_log10_p)*0.98, 
           label = paste0("n=", n_down), color = COL_HEALTHY, size = 3.2, fontface = "bold") +
  annotate("text", x = max(df$z_score)*0.8, y = max(df$neg_log10_p)*0.98, 
           label = paste0("n=", n_up), color = COL_DCM, size = 3.2, fontface = "bold") +
  
  # Etiquetas de quinasas (ajustadas para el tamaño de punto nuevo)
  geom_text_repel(data = top_labels, aes(label = kinase),
                  size = 2.1, color = "black", fontface = "italic",
                  box.padding = 0.3, point.padding = 0.4,
                  segment.linewidth = 0.2, min.segment.length = 0) +
  
  # Estética y Tema limpio
  labs(x = "Kinase Activity (z-score)", y = "-log10(p-value)") +
  theme_void() + 
  theme(
    axis.title       = element_text(size = 7.5, color = "grey20"),
    axis.text        = element_text(size = 6.5, color = "grey30"),
    axis.line        = element_line(color = "black", linewidth = 0.35),
    axis.ticks       = element_line(color = "black", linewidth = 0.25),
    plot.margin      = margin(8, 8, 8, 8),
    legend.position  = "none"
  )

# ==============================================================================
# 3. EXPORTACIÓN
# ==============================================================================
ggsave("KSEA_Volcano_LargePoints.pdf", p, width = w_in, height = h_in, units = "in", device = "pdf")
ggsave("KSEA_Volcano_LargePoints.jpg", p, width = w_in, height = h_in, units = "in", dpi = 300)





#####################
# GSEA HALLMARKS KINASE ENRICHMENT
####################
library(ggplot2)
library(dplyr)
library(stringr)

# ==============================================================================
# 1. CONFIGURACIÓN Y DATOS
# ==============================================================================
COL_DCM     <- "#F5CD6A"  # Oro (Disease)
COL_HEALTHY <- "#3A2044"  # Morado (Healthy)

# Dimensiones exactas (169.05 x 174.63 pt)
w_in <- 169.0586 / 72
h_in <- 174.6366 / 72

df <- read.csv("fgsea_hallmarks_results.csv")

# Preparar datos: Filtro por pval < 0.05
df_plot <- df %>%
  filter(pval < 0.05) %>%
  mutate(
    pathway_clean = str_remove(pathway, "HALLMARK_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " "),
    pathway_clean = str_to_title(pathway_clean),
    log_p = -log10(pval)
  ) %>%
  arrange(NES)

df_plot$pathway_clean <- factor(df_plot$pathway_clean, levels = df_plot$pathway_clean)

# ==============================================================================
# 2. GENERAR LOLLIPOP PLOT CON LEYENDAS EXPLICATIVAS
# ==============================================================================
p <- ggplot(df_plot, aes(x = NES, y = pathway_clean)) +
  # Línea del lollipop
  geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean, 
                   color = direction), linewidth = 0.6) +
  
  # Círculo (Tamaño = p-value)
  geom_point(aes(color = direction, size = log_p), stroke = 0) +
  
  # Estrellas de significancia
  geom_text(aes(x = if_else(NES > 0, NES + 0.25, NES - 0.25), label = sig_label), 
            size = 2, vjust = 0.8) +
  
  # Configuración de Colores (Leyenda de Condición)
  scale_color_manual(values = c("Enriched in Disease" = COL_DCM, 
                                "Enriched in Healthy" = COL_HEALTHY),
                     name = "Condition") +
  
  # Configuración de Tamaños (Leyenda de p-value)
  scale_size_continuous(range = c(1.5, 4.5), 
                        name = expression("-log"[10]*"(p)")) +
  
  # Línea central
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
  # Estética y Leyendas
  labs(x = "NES", y = NULL) +
  theme_minimal() +
  theme(
    # Quitar rejillas y fondo
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    # Ejes
    axis.title.x     = element_text(size = 5.5, face = "bold"),
    axis.text.x      = element_text(size = 5),
    axis.text.y      = element_text(size = 4, color = "black", hjust = 1),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    
    # LEYENDAS (Ajustadas para tamaño pequeño)
    legend.title     = element_text(size = 4.5, face = "bold"),
    legend.text      = element_text(size = 4),
    legend.key.size  = unit(0.2, "cm"),
    legend.position  = "right",
    legend.margin    = margin(0, 0, 0, 0),
    
    plot.margin      = margin(5, 2, 5, 2)
  )

# ==============================================================================
# 3. EXPORTACIÓN
# ==============================================================================
ggsave("GSEA_Lollipop_Final_Legend.pdf", p, width = w_in, height = h_in, units = "in", device = "pdf")
ggsave("GSEA_Lollipop_Final_Legend.jpg", p, width = w_in, height = h_in, units = "in", dpi = 300)



#######################
## reactome
#####################

library(ggplot2)
library(dplyr)
library(stringr)

# ==============================================================================
# 1. CONFIGURACIÓN Y DATOS
# ==============================================================================
COL_DCM     <- "#F5CD6A"  # Oro (Disease)
COL_HEALTHY <- "#3A2044"  # Morado (Healthy)

# Dimensiones exactas (169.0586 x 174.6366 pt)
w_in <- 169.0586 / 72
h_in <- 174.6366 / 72

# Lista de términos solicitados
target_pathways <- c(
  "REACTOME_PI3K_EVENTS_IN_ERBB4_SIGNALING",
  "REACTOME_MAPK1_MAPK3_SIGNALING",
  "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
  "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
  "REACTOME_MAPK3_ERK1_ACTIVATION",
  "REACTOME_PI_3K_CASCADE_FGFR2",
  "REACTOME_ERK_MAPK_TARGETS",
  "REACTOME_PYRUVATE_METABOLISM",
  "REACTOME_SIGNALING_BY_ACTIVIN",
  "REACTOME_SIGNALLING_TO_ERKS"
)

df <- read.csv("fgsea_reactome_results.csv")

# Preparar datos
df_plot <- df %>%
  filter(pathway %in% target_pathways) %>%
  mutate(
    # Limpiar nombres: quitar "REACTOME_" y poner en formato título
    pathway_clean = str_remove(pathway, "REACTOME_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " "),
    pathway_clean = str_to_title(pathway_clean),
    # Significancia para el tamaño (-log10 p-value)
    log_p = -log10(pval)
  ) %>%
  arrange(NES)

# Ordenar el eje Y por el valor de NES
df_plot$pathway_clean <- factor(df_plot$pathway_clean, levels = df_plot$pathway_clean)

# ==============================================================================
# 2. GENERAR LOLLIPOP PLOT
# ==============================================================================
p <- ggplot(df_plot, aes(x = NES, y = pathway_clean)) +
  # Línea del lollipop (brazo)
  geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean, 
                   color = direction), linewidth = 0.6) +
  
  # Círculo del lollipop (Tamaño = log_p, Color = direction)
  geom_point(aes(color = direction, size = log_p), stroke = 0) +
  
  # Colores manuales
  scale_color_manual(values = c("Enriched in Disease" = COL_DCM, 
                                "Enriched in Healthy" = COL_HEALTHY),
                     name = "Condition") +
  
  # Escala de tamaño para el p-value
  scale_size_continuous(range = c(1.5, 4.5), 
                        name = expression("-log"[10]*"(p)")) +
  
  # Línea central de referencia (NES = 0)
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
  # Etiquetas y Tema Limpio (Sin rejillas)
  labs(x = "NES", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    # Ejes
    axis.title.x     = element_text(size = 5.5, face = "bold"),
    axis.text.x      = element_text(size = 5),
    axis.text.y      = element_text(size = 4, color = "black", hjust = 1),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    
    # LEYENDAS (Ajustadas para el tamaño reducido)
    legend.title     = element_text(size = 4.5, face = "bold"),
    legend.text      = element_text(size = 4),
    legend.key.size  = unit(0.2, "cm"),
    legend.position  = "right",
    legend.margin    = margin(0, 0, 0, 0),
    
    plot.margin      = margin(5, 2, 5, 2)
  )

# ==============================================================================
# 3. EXPORTACIÓN
# ==============================================================================
ggsave("Reactome_Lollipop_Pvalue.pdf", p, width = w_in, height = h_in, units = "in", device = "pdf")
ggsave("Reactome_Lollipop_Pvalue.jpg", p, width = w_in, height = h_in, units = "in", dpi = 300)






###########
## HEATMAP
##########
# 1. Cargar librerías
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# 2. Leer datos
df <- read.csv("kinase_activity_all_methods.csv")

# 3. Filtrar y preparar la matriz
# Usamos p_ksea < 0.05 y seleccionamos las columnas de actividad
df_sig <- df %>%
  filter(p_ksea < 0.05) %>%
  select(kinase, activity_ksea, activity_ulm, activity_wmean) %>%
  rename(KSEA = activity_ksea, ULM = activity_ulm, WMEAN = activity_wmean)

# Convertir a matriz y poner nombres de las quinasas en las filas
mat <- as.matrix(df_sig[,-1])
rownames(mat) <- df_sig$kinase

# 4. Configuración de colores (Oro para +, Morado para -)
COL_DCM     <- "#F5CD6A"
COL_HEALTHY <- "#3A2044"

# Definimos la escala: Morado (mínimo) -> Blanco (0) -> Oro (máximo)
limit <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-limit, 0, limit), c(COL_HEALTHY, "white", COL_DCM))

# 5. Crear el Heatmap
ht <- Heatmap(mat, 
              name = "Activity KSEA\n(Z-score)", 
              col = col_fun,
              cluster_columns = FALSE,         # Mantener orden de los métodos
              show_row_dend = TRUE,            # Agrupar quinasas similares
              row_names_gp = gpar(fontsize = 6), # Fuente pequeña para que quepan todas
              column_names_gp = gpar(fontsize = 9, fontface = "bold"),
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 0.5), # Línea blanca entre celdas
              column_title = "Kinase Activity (p_ksea < 0.05)")

# 6. GUARDAR EL ARCHIVO (Corregido para evitar el error de carpeta)
# Guardamos directamente en el directorio actual para evitar errores de ruta
png("Kinase_Activity_Heatmap.png", width = 800, height = 1600, res = 150)
draw(ht)
dev.off()

pdf("Kinase_Activity_Heatmap.pdf", width = 5, height = 12)
draw(ht)
dev.off()

cat("Heatmap guardado exitosamente como 'Kinase_Activity_Heatmap.png' y '.pdf'\n")


######
## heaetmap compact
####

# 1. Cargar librerías
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(circlize)

# 2. Cargar y filtrar datos
df <- read.csv("kinase_activity_all_methods.csv")

# Filtramos por p_ksea < 0.05 y seleccionamos las columnas de comparación
df_sig <- df %>%
  filter(p_ksea < 0.05) %>%
  select(kinase, activity_ksea, activity_ulm, activity_wmean) %>%
  rename(KSEA = activity_ksea, ULM = activity_ulm, WMEAN = activity_wmean)

# Convertir a matriz
mat <- as.matrix(df_sig[,-1])
rownames(mat) <- df_sig$kinase

# 3. Configurar escala de color MAGMA
# Mapeamos los valores al rango completo de la paleta Magma
col_fun <- colorRamp2(
  seq(min(mat, na.rm=TRUE), max(mat, na.rm=TRUE), length.out = 100), 
  viridis::magma(100)
)

# 4. Crear el Heatmap Ampliado
# Restauramos el dendrograma y aumentamos el tamaño de letra
ht <- Heatmap(mat, 
              name = "Activity (Z-score)", 
              col = col_fun,
              cluster_columns = FALSE,         # Mantener orden: KSEA, ULM, WMEAN
              show_row_dend = TRUE,            # Mostrar agrupamiento de quinasas
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 8), # Letra legible
              column_names_gp = gpar(fontsize = 10, fontface = "bold"),
              column_names_rot = 0,            # Nombres de columnas horizontales
              
              # Ajuste de celdas
              column_title = "Kinase Activity Comparison (p_ksea < 0.05)",
              border = TRUE,
              rect_gp = gpar(col = "white", lwd = 0.5))

# 5. EXPORTACIÓN EN FORMATO GRANDE
# Dimensiones recomendadas para que no se vea "feote" ni amontonado
# 7 pulgadas de ancho x 14 de alto para acomodar la lista larga de quinasas
pdf("Kinase_Comparison_Heatmap_Magma_Large.pdf", width = 7, height = 14)
draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

jpeg("Kinase_Comparison_Heatmap_Magma_Large.jpg", width = 7, height = 14, units = "in", res = 300)
draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

cat("Heatmap ampliado generado exitosamente.\n")
cat("Se han incluido las comparaciones de KSEA, ULM y WMEAN para quinasas con p < 0.05.\n")
