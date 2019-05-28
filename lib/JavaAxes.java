import java.awt.*;
import java.awt.geom.Line2D;
import javax.swing.JLabel;
import javax.swing.JPanel;

// this class is based on the code written in https://stackoverflow.com/a/18413639
public class JavaAxes extends JPanel {

    private int padding = 20;
    private int labelPadding = 25;
	private Color lineColor = Color.BLUE;
    private Color gridColor = new Color(200, 200, 200, 200);
    private int tickLength = 4;	
	private int xticksNumber = 4;
    private int yticksNumber = 6;
    private double[] xdata, ydata;
	private String xlabel, ylabel;
	private double xfactor = 1;
	private double yfactor = 1;
	private double xmin, xmax, ymin, ymax;
	
	public JavaAxes() {
		setXLim(new double[]{0, 1});
		setYLim(new double[]{0, 1});
		setAxisLabels("", "");
	}

    public JavaAxes(double[] x,double[] y) {        
		xdata = x;
		ydata = y;
		
		// define the axis limits according the data
		setXLim(new double[0]);
		setYLim(new double[0]);
		
		setAxisLabels("", "");
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // draw white background
		int panelHeight = getHeight();
		int panelWidth = getWidth();
        g2.setColor(Color.WHITE);
        g2.fillRect(padding + labelPadding, padding, panelWidth - 2 * padding - labelPadding, panelHeight - 2 * padding - labelPadding);        
		
		FontMetrics metrics = g2.getFontMetrics();		
		
		// get graphic scale
		double xScale = ((double) panelWidth - 2 * padding - labelPadding) / (xmax - xmin);
		double yScale = ((double) panelHeight - 2 * padding - labelPadding) / (ymax - ymin);		

        // draw tick marks and grid lines for y-axis
		String yticklabel;
		int ytickLabelWidth;
		double dy = (ymax - ymin) * yScale / yticksNumber;
		double y0, y1;		
		double x0 = padding + labelPadding;
		double x1 = tickLength + padding + labelPadding;		
        for (int i = 0; i < yticksNumber + 1; i++) {
            y0 = panelHeight - (i * dy + padding + labelPadding);
            y1 = y0;
			
			g2.setColor(gridColor);
			g2.draw(new Line2D.Double(x0 + 1 + tickLength, y0, panelWidth - padding, y1));
			
			g2.setColor(Color.BLACK);
			yticklabel = (int) Math.round((ymin + i * (ymax - ymin) / yticksNumber) * (100 / yfactor)) / 100.0 + "";
			ytickLabelWidth = metrics.stringWidth(yticklabel);
			g2.drawString(yticklabel, (int) x0 - ytickLabelWidth - 5, (int) y0 + (metrics.getHeight() / 2) - 3);

			g2.draw(new Line2D.Double(x0, y0, x1, y1));
        }

		// draw tick marks and grid lines for x-axis
		String xticklabel;
		int xtickLabelWidth;
		double dx = (xmax - xmin) * xScale / xticksNumber;		
		y0 = panelHeight - padding - labelPadding;
		y1 = y0 - tickLength;
		for (int i = 0; i < xticksNumber + 1; i++) {			
			x0 = padding + labelPadding + i * dx;
			x1 = x0;
			
			g2.setColor(gridColor);
			g2.draw(new Line2D.Double(x0, y0 - 1 - tickLength, x1, padding));
			
			g2.setColor(Color.BLACK);
			xticklabel = (int) Math.round((xmin + i * (xmax - xmin) / xticksNumber) * (100 / xfactor)) / 100.0 + "";
			xtickLabelWidth = metrics.stringWidth(xticklabel);
			g2.drawString(xticklabel, (int) x0 - xtickLabelWidth / 2, (int) y0 + metrics.getHeight() + 3);

			g2.draw(new Line2D.Double(x0, y0, x1, y1));			
		}

        // draw axes
        g2.drawLine(padding + labelPadding, (int) y0, padding + labelPadding, padding);
        g2.drawLine(padding + labelPadding, (int) y0, panelWidth - padding, (int) y0);
        
		// draw line plot of the data
        g2.setColor(lineColor);
		if (xdata != null) {			
			int n = xdata.length;
			double[] xp = new double[n];
			double[] yp = new double[n];			

			for (int i = 0; i < n; i++) {
				xp[i] =  (xdata[i] - xmin) * xScale + padding + labelPadding;
				yp[i] =  (ymax - ydata[i]) * yScale + padding;
			}
		
			for (int i = 0; i < n - 1; i++)
				g2.draw(new Line2D.Double(xp[i], yp[i], xp[i + 1], yp[i + 1]));
		}
    }
	
    private double getMin(double[] vector) {
        double min = vector[0];
        for (int i = 1; i < vector.length; i++) {
            min = Math.min(min, vector[i]);
        }
        return min;
    }

    private double getMax(double[] vector) {
        double max = vector[0];
        for (int i = 1; i < vector.length; i++) {
            max = Math.max(max, vector[i]);
        }
        return max;
    }
	
	private void defineXFactor() {		
		double powerX = Math.log10(Math.max(Math.abs(xmin), Math.abs(xmax)));
		if (powerX < -2 || powerX > 3)
			xfactor = Math.pow(10, Math.signum(powerX) * Math.ceil(Math.abs(powerX)));
		else		
			xfactor = 1;
	}
	
	private void defineYFactor() {
		double powerY = Math.log10(Math.max(Math.abs(ymin), Math.abs(ymax)));
		if (powerY < -2 || powerY > 3)			
			yfactor = Math.pow(10, Math.signum(powerY) * Math.ceil(Math.abs(powerY)));
		else		
			yfactor = 1;
	}

    public void plot(double[] x, double[] y) {
		xdata = x;
		ydata = y;
		setXLim(new double[0]);
		setYLim(new double[0]);
		setAxisLabels(xlabel, ylabel);
    }	
	
	public void setXLim(double[] xlim) {
		if (xlim.length == 0) { // define the x-axis limits according the data
			xmin = getMin(xdata);
			xmax = getMax(xdata);
			if (xmin == xmax) {
				if (xmin != 0) {
					xmin = 0;
					xmax = 2 * xmax;
				}
				else {
					xmin = -0.5;
					xmax = 0.5;
				}
			}
			defineXFactor();
		}
		else {
			xmin = xlim[0];
			xmax = xlim[1];
			defineXFactor();
			setAxisLabels(xlabel, ylabel);
		}		
	}
	
	public void setYLim(double[] ylim) {
		if (ylim.length == 0) { // define the y-axis limits according the data
			ymin = getMin(ydata);
			ymax = getMax(ydata);
			if (ymin == ymax) {
				if (ymin != 0) {
					ymin = 0;
					ymax = 2 * ymax;
				}
				else {
					ymin = -0.5;
					ymax = 0.5;
				}
			}
			defineYFactor();
		}
		else {
			ymin = ylim[0];
			ymax = ylim[1];
			defineYFactor();
			setAxisLabels(xlabel, ylabel);
		}
	}	
	
	public void setAxisLabels(String x_label, String y_label) {				
		// remove old labels (and everything else...)
		removeAll();
		
		xlabel = x_label;
		ylabel = y_label;
		
		// create layout to organize the labels
		setLayout(new GridBagLayout());		
		GridBagConstraints grid = new GridBagConstraints();		
		
		// include a y-axis factor if necessary
		int flag_yfactor = 0;		
		JLabel yFactorLabel = new JLabel("");
		if (yfactor != 1) {
			flag_yfactor = 1;
			grid.insets = new Insets(0, padding + labelPadding, 0, 0);
			grid.anchor = GridBagConstraints.NORTHWEST;
			grid.weightx = 1;			
			yFactorLabel.setText("<html>&times;10<sup>" + Math.round(Math.log10(yfactor)));
			add(yFactorLabel, grid);
			grid.gridy = 1;
		}
		
		// add the y-axis label
		grid.insets = new Insets(padding, 0, padding + labelPadding + flag_yfactor * 22, 0);
		grid.anchor = GridBagConstraints.WEST;
		grid.weightx = 1;	
		JLabel YLabel = new JLabel(ylabel);
		add(YLabel, grid);
		
		// include a x-axis factor if necessary
		int flag_xfactor = 0;
		JLabel xFactorLabel = new JLabel("");
		if (xfactor != 1) {
			flag_xfactor = 1;
			xFactorLabel.setText("<html>&times;10<sup>" + Math.round(Math.log10(xfactor)));
		}	
		
		// add the x-axis label
		int delta = labelPadding - YLabel.getPreferredSize().width;
		grid.weighty = 1;
		grid.insets = new Insets(0, padding + delta + flag_xfactor * xFactorLabel.getPreferredSize().width - flag_yfactor * (yFactorLabel.getPreferredSize().width + padding + delta), padding + labelPadding - 40, 0);
		grid.anchor = GridBagConstraints.SOUTH;
		add(new JLabel(xlabel), grid);		
		
		grid.anchor = GridBagConstraints.SOUTHEAST;
		grid.insets = new Insets(0, 0, padding + labelPadding - 40, padding);
		add(xFactorLabel, grid);		
		
		revalidate();
        repaint();		
	}

	public void setPadding(int padd) {
		padding = padd;
        setAxisLabels(xlabel, ylabel);        
    }
	
	public void setLabelPadding(int padd) {
		labelPadding = padd;
        setAxisLabels(xlabel, ylabel);
    }
	
	public void setXTicksNumber(int num) {
		xticksNumber = num;
		revalidate();
        repaint();
	}
	
	public void setYTicksNumber(int num) {
		yticksNumber = num;
		revalidate();
        repaint();
	}
}
