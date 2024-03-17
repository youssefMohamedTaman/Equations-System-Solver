module FX_Project2 {
	requires javafx.controls;
	requires javafx.fxml;
	requires commons.math3;
	
	opens application to javafx.graphics, javafx.fxml;
}
