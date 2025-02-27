##http://www.cnblogs.com/xudongliang/p/7884667.html-original code for flower plot
##https://www.biostars.org/p/144778/

flower_plot <- function(sample, value, start, a, b,  
                        ellipse_col = rgb(135, 206, 235, 150, max = 255), 
                        circle_col = rgb(0, 162, 214, max = 255),
                        circle_text_cex = 1, labels=labels) {
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type="n")
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                          col = ellipse_col,
                          border = ellipse_col,
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t]
    )
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = circle_text_cex
      )
      
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = circle_text_cex
      )
    }           
  })
  plotrix::draw.circle(x = 5, y = 5, r = 0.65, col = circle_col, border = circle_col)
  text(x = 5, y = 5, labels=labels)
}


flower_plot2 <- function(sample, value, start, a, b,  
                         ellipse_col = rgb(135, 206, 235, 150, max = 255), 
                         circle_col = rgb(0, 162, 214, max = 255),
                         circle_text_cex = 1, labels=labels) {
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type="n")
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    ellipse_col <- ellipse_col[t]
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                          col = ellipse_col,
                          border = ellipse_col,
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t]
    )
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = circle_text_cex
      )
      
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = circle_text_cex
      )
    }           
  })
  plotrix::draw.circle(x = 5, y = 5, r = 0.45, col = circle_col, border = circle_col )
  
  # tune location by x and y.
  text(x = 5, y = 5, labels=labels)
}



flower_plot2(c(
                 'S. Anatum',
                'S. Berta',
              'S. Braenderup',
                 'S. Cerro','S. Enteritidis',
                'S. Give',
                'S. Hartford',
               'S. Montevideo',
               'S. Muenchen',
                'S. Newport',
                'S. Oranienburg'), 
              c( 117, 43, 304, 746, 153, 472, 626, 311, 110, 946, 510), 90, 0.45, 2, labels="3436",
              ellipse_col = rainbow(15, alpha = 0.5),
             circle_col = "beige")
            
    
