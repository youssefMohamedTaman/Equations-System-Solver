package application;

import SystemSolver.Solver;
import SystemSolver.secant;
import SystemSolver.Bisection;
import SystemSolver.FalsePosition;
import SystemSolver.FixedPoint;
import SystemSolver.GaussElimination;
import SystemSolver.GaussSeidel;
import SystemSolver.GuassJordanElimination;
import SystemSolver.Jacobi;
import SystemSolver.LUCholesky;
import SystemSolver.LUCrout;
import SystemSolver.LUDoolittle;
import SystemSolver.ModifiedNewtonRaphson;
import SystemSolver.NewtonRaphason;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.geometry.Pos;
import javafx.scene.control.TextField;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.VBox;

import java.net.URL;
import java.util.ResourceBundle;
import java.util.function.Function;

public class MainSceneController implements Initializable {
	int N = 0;
	int significantDigits = 10;
	String method = "";
	int iter = 100;
	double er = 0;
	@FXML
	private ComboBox<String> comBox;
	@FXML
	private TextField tfNumOfEq;
	@FXML
	private TextField signDigits;
	@FXML
	private AnchorPane ap;
    @FXML
    private ScrollPane sp;
    @FXML
	private VBox vb;
	@Override 
	public void initialize(URL url,ResourceBundle resourseBundle) {
		comBox.setItems(FXCollections.observableArrayList
				("Gauss Elimination","Gauss-Jordan",
				"LU Decomposition(Crout)","LU Decomposition(Cholesky)",
				"LU Decomposition(Doolittle)","Gauss-Seidel",
				"Jacobi-Iteration","Bisection","False-Position","Fixed Point","Secant",
				"Newton Raphason","Modified Newton Raphason 1","Modified Newton Raphason 2"));
	}
	@FXML
	public void btnOKClicked(ActionEvent event) {
		// reading significant figures input
		if(!signDigits.getText().equals(""))
			significantDigits = Integer.parseInt(signDigits.getText());
		else 
			significantDigits = 10; 	// if there is no input initializing it to 5
		//creating the coefficient matrix
		String Num = tfNumOfEq.getText();
		if(Num!=null && !Num.equals("") && !Num.equals("0")){
			//if(N != Integer.parseInt(Num)){
			
			N = Integer.parseInt(Num);
			if(method.equals("JI") || method.equals("GS"))
				create2(N);
			
			if(method.equals("FP") || method.equals("BIS") ||
					   method.equals("FXP") ||method.equals("SEC")||
					   method.equals("NR")||method.equals("MNR")||method.equals("MNR2")) {
						create3(N);
						}
			else create(N);
			}
		//}	
	}
	// creating the primary objects in the UI
	public void create(int n) {
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if((node.getId().equals("=")) || 
				(node.getId().equals("solve")) || 
				(node.getId().charAt(0)=='a') || 
				(node.getId().charAt(0)=='x') || 
				(node.getId().charAt(0)=='b') ||
				(node.getId().equals("rt")) ||
				(node.getId().equals("Lbound")) ||
				(node.getId().equals("Ubound")) ||
				(node.getId().equals("BL")) || 
				node.getId().contains("dl") ||
				node.getId().contains("cc") ||
				node.getId().equals("mm") ||
				node.getId().equals("reachedIt")
				)
				return true;
			return false;
		});
		if(method.equals("FP") || method.equals("BIS") ||  method.equals("FXP")) {
			return;
		}
	    int tfWidth = 30;
	    int tfHeight = 28;
	    int startX = 35;
	    int startY = 200;

	    for (int i = 0; i < n; i++) {
	        int x = startX;
	        int y = startY + i * (tfHeight + 10); 	// Adjust the vertical position
            // creating the coefficients matrix(n*n) A 
	        for (int j = 0; j < n; j++) {
	            String w = Integer.toString(i);
	            String q = Integer.toString(j);
	            String tfId = "a" + w + q;
	            String labelId = "x" + w + q;
	            TextField TF = new TextField();
	            Label L = new Label();
	            TF.setId(tfId);
	            TF.setLayoutX(x);
	            TF.setLayoutY(y);
	            TF.setPrefHeight(tfHeight);
	            TF.setPrefWidth(tfWidth);
	            TF.setStyle("-fx-background-color:  #33887799;");
	            x += tfWidth + 10;
	            L.setId(labelId);
	            L.setText("X"+q);
	            L.setLayoutX(x);
	            L.setLayoutY(y);
	            L.setPrefHeight(tfHeight);
	            L.setPrefWidth(tfWidth);
	            x += tfWidth + 10; 	// Adjust the horizontal position
	            ap.getChildren().addAll(TF, L);
	        }
	        x -= 10; 	// Adjust the horizontal position
	        
	        // creating the constants vector(n) B 
	        TextField TF = new TextField();
            Label L = new Label();
            L.setId("=");
            L.setText("=");
            L.setLayoutX(x);
            L.setLayoutY(y);
            L.setPrefHeight(tfHeight);
            L.setPrefWidth(tfWidth);
            x += 20; 	// Adjust the horizontal position
            TF.setId("b"+Integer.toString(i));
            TF.setLayoutX(x);
            TF.setLayoutY(y);
            TF.setPrefHeight(tfHeight);
            TF.setPrefWidth(tfWidth);
            TF.setStyle("-fx-background-color:  #33887799;");
            ap.getChildren().addAll(TF, L);   
	    }
	    // creating a button to solve the system we created
	    Button b = new Button();
	    b.setId("solve");
	    b.setText("Solve");
	    b.setLayoutX(2*n*40 + 35 );
	    b.setLayoutY(n*28 + 200 + 10*n);
	    b.setStyle("-fx-font-weight: bold;");
	    b.setOnAction(event -> solve());
	    b.hoverProperty().addListener(l->{
	    	b.setStyle("-fx-background-color:  #33887799; -fx-font-weight: bold;");
	    });

	    b.setOnMouseMoved(m->{
	    b.setStyle("-fx-background-color: #3874; -fx-font-weight: bold;");
	    });
	    ap.getChildren().add(b);
	    vb.setPrefSize(80*n+100, 40*n+300);
	    
	}
	//creating the secondary objects in the UI which will be created only while using Jacobi and gauss_siedel 
	public void create2(int n) {
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if(node.getId().contains("ix") || 
					node.getId().contains("ia")||
					node.getId().equals("numOfIt")||
					node.getId().equals("stopCondition")||
					node.getId().equals("error")||
					node.getId().equals("v1")||
					node.getId().equals("v2")
					)
				return true;
			return false;
		});
		Label La = new Label();
		La.setId("ix");
        La.setText("Initial guess:");
        La.setLayoutX(16);
        La.setLayoutY(60);
        La.setPrefHeight(28);
        La.setPrefWidth(80);
        La.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
        La.setAlignment(Pos.CENTER);
        ap.getChildren().add(La);
		int x = 95 ;
		// vreating intial guess vector
		for(int i=0;i<n;i++) {
			String w = Integer.toString(i);
			TextField TF = new TextField();
            Label L = new Label();
            L.setId("ix"+w);
            L.setText("X"+w+" = ");
            L.setLayoutX(x);
            L.setLayoutY(60);
            L.setPrefHeight(28);
            L.setPrefWidth(35);
            x +=35;
            TF.setId("ia"+w);
            TF.setLayoutX(x);
            TF.setLayoutY(60);
            TF.setPrefHeight(28);
            TF.setPrefWidth(35);
            TF.setStyle("-fx-background-color:  #33887799;");
            x += 45;
            
            ap.getChildren().addAll(TF,L);
		}
		// creating maximum number of iteration text input 
		Label L1 = new Label();
		Label L2 = new Label();
		Label L3 = new Label();
		TextField T1 = new TextField();
		TextField T2 = new TextField();
		L1.setId("stopCondition");
        L1.setText("Stoping Condition:");
        L1.setLayoutX(0);
        L1.setLayoutY(90);
        L1.setPrefHeight(28);
        L1.setPrefWidth(150);
        L1.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
        L1.setAlignment(Pos.CENTER);
        
        L2.setId("numOfIt");
        L2.setText("Number Of Iterations:");
        L2.setLayoutX(8);
        L2.setLayoutY(120);
        L2.setPrefHeight(28);
        L2.setPrefWidth(150);
        L2.setStyle(" -fx-font-weight: bold; -fx-background-radius:5;");
        L2.setAlignment(Pos.CENTER);
        
        T1.setId("v1");
        T1.setLayoutX(154);
        T1.setLayoutY(120);
        T1.setPrefHeight(28);
        T1.setPrefWidth(45);
        T1.setStyle("-fx-background-color:  #33887799;");
        //creating relative error text input 
        L3.setId("error");
        L3.setText("Absolute Relative Error:");
        L3.setLayoutX(214);
        L3.setLayoutY(120);
        L3.setPrefHeight(28);
        L3.setPrefWidth(200);
        L3.setStyle(" -fx-font-weight: bold; -fx-background-radius:5;");
        L3.setAlignment(Pos.CENTER);
        
        T2.setId("v2");
        T2.setLayoutX(389);
        T2.setLayoutY(120);
        T2.setPrefHeight(28);
        T2.setPrefWidth(60);
        T2.setStyle("-fx-background-color: #33887799;");
        
        ap.getChildren().addAll(L1,L2,L3,T1,T2);
		
	}
	//creating objects of non linear equations
	public void create3(int n){
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if((node.getId().equals("=")) || 
				(node.getId().equals("solve")) || 
				(node.getId().charAt(0)=='a') || 
				(node.getId().charAt(0)=='x') || 
				(node.getId().charAt(0)=='b') ||
				(node.getId().equals("rt")) ||
				(node.getId().equals("Lbound")) ||
				(node.getId().equals("Ubound")) ||
				(node.getId().equals("BL")) ||
				node.getId().equals("numOfIt")||
				node.getId().equals("stopCondition")||
				node.getId().equals("error")||
				node.getId().equals("v1")||
				node.getId().equals("v2") ||
				node.getId().contains("dl") ||
				node.getId().contains("cc") ||
				node.getId().equals("mm") ||
				node.getId().equals("reachedIt")
				)
				return true;
			return false;
		});
	    int tfWidth = 60;
	    int tfHeight = 28;
	    int startX = 35;
	    
	         create4();
	         int t = startX;
	         int len = t;
	         if(method.equals("FXP") || method.equals("NR")|| method.equals("MNR") || method.equals("MNR2")) {
		            TextField TF1 = new TextField();
		            Label L4 = new Label();
		            L4.setId("BL");
		            L4.setText("Intial guess X0:");
		            L4.setLayoutX(startX);
		            L4.setLayoutY(150);
		            L4.setPrefHeight(tfHeight);
		            L4.setPrefWidth(90);
		            L4.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
		            TF1.setId("Lbound");
		            TF1.setLayoutX(startX+100);
		            TF1.setLayoutY(150);
		            TF1.setPrefHeight(tfHeight);
		            TF1.setPrefWidth(tfWidth);
		            TF1.setStyle("-fx-background-color:  #33887799;");
		            ap.getChildren().addAll(L4,TF1);
		            if(method.equals("MNR2")) {
		            	  TextField TF2 = new TextField();
				            Label L5 = new Label();
				            L5.setId("BL");
				            L5.setText("multiplicity(m) =");
				            L5.setLayoutX(startX+200);
				            L5.setLayoutY(150);
				            L5.setPrefHeight(tfHeight);
				            L5.setPrefWidth(120);
				            L5.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
				            TF2.setId("mm");
				            TF2.setLayoutX(startX+310);
				            TF2.setLayoutY(150);
				            TF2.setPrefHeight(tfHeight);
				            TF2.setPrefWidth(tfWidth);
				            TF2.setStyle("-fx-background-color:  #33887799;");
				            ap.getChildren().addAll(L5,TF2);
		            }
	         }
		     else if(method.equals("BIS") || method.equals("FP") || method.equals("SEC")){
				    	for (int i = 1; i >=0; i--) {
				    		String tfId;
				    		if(i==1)
				    			tfId = "Lbound";
				    		else
				    			tfId = "Ubound";
				            String labelId = "BL" ;
				            TextField TF = new TextField();
				            Label L = new Label();
				            L.setId(labelId);
				            if(i==1)
				            	L.setText("Lower guess Xl:");
				            else
				            	L.setText("Upper guess Xu:");
				            L.setLayoutX(len);
				            L.setLayoutY(150);
				            L.setPrefHeight(tfHeight);
				            L.setPrefWidth(95);
				            L.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
				            len += 100;
				            TF.setId(tfId);
				            TF.setLayoutX(len);
				            TF.setLayoutY(150);
				            TF.setPrefHeight(tfHeight);
				            TF.setPrefWidth(tfWidth);
				            TF.setStyle("-fx-background-color:  #33887799;");
				            len += 100;
				     
				            ap.getChildren().addAll(TF, L);
					    	}
		     }
	        	 for(int i=0;i<n;i++) {
	        		 String w = Integer.toString(i);
			            Label L = new Label();
			            L.setId("x"+w);
			            L.setText("coeff:");
			            L.setLayoutX(t);
			            L.setLayoutY(190);
			            L.setPrefHeight(tfHeight);
			            L.setPrefWidth(30);
			            t += 40;
			            TextField TF = new TextField();
			            TF.setId("a"+w);
			            TF.setLayoutX(t);
			            TF.setLayoutY(190);
			            TF.setPrefHeight(tfHeight);
			            TF.setPrefWidth(tfWidth);
			            TF.setStyle("-fx-background-color:  #33887799;");
			            t += tfWidth + 10;
	        		 ComboBox<String> dl = new ComboBox<>();
	        		 dl.setId("dl"+w);
	        		 dl.setPromptText("Select function");
	        		 dl.getItems().addAll("x^c", "e^(cx)", "sin(cx)", "cos(cx)");
	        		 dl.setLayoutX(t);
	        		 dl.setLayoutY(190);
	        		 dl.setPrefHeight(29);
	        		 dl.setPrefWidth(140);
	        		 t += 150;
	        		 	Label LL = new Label();
			            LL.setId("x"+w);
			            LL.setText("c=");
			            LL.setLayoutX(t);
			            LL.setLayoutY(190);
			            LL.setPrefHeight(tfHeight);
			            LL.setPrefWidth(15);
			            t += 25;
			            TextField TFF = new TextField();
			            TFF.setId("cc"+w);
			            TFF.setLayoutX(t);
			            TFF.setLayoutY(190);
			            TFF.setPrefHeight(tfHeight);
			            TFF.setPrefWidth(tfWidth);
			            TFF.setStyle("-fx-background-color:  #33887799;");
			            t += tfWidth + 10;
			            ap.getChildren().addAll(dl,L,TF,LL,TFF);
	        	 }

	    Button b = new Button();
	    b.setId("solve");
	    b.setText("Solve");
	    b.setLayoutX(startX);
	    b.setLayoutY(240);
	    b.setStyle("-fx-font-weight: bold;");
	    b.setOnAction(event -> solve2(n));
	    b.hoverProperty().addListener(l->{
	    	b.setStyle("-fx-background-color:  #33887799; -fx-font-weight: bold;");
	    });
	    b.setOnMouseMoved(m->{
	    b.setStyle("-fx-background-color: #3874; -fx-font-weight: bold;");
	    });
	    ap.getChildren().add(b);
	    vb.setPrefWidth(350*n+20*n+200);
	}
	public void create4() {
 		Label L1 = new Label();
 		Label L2 = new Label();
 		Label L3 = new Label();
 		TextField T1 = new TextField();
 		TextField T2 = new TextField();
 		
 		 L1.setId("stopCondition");
         L1.setText("Stoping Condition:");
         L1.setLayoutX(0);
         L1.setLayoutY(70);
         L1.setPrefHeight(28);
         L1.setPrefWidth(150);
         L1.setStyle("-fx-font-weight: bold; -fx-background-radius:5;");
         L1.setAlignment(Pos.CENTER);
         // creating maximum number of iteration text input
         L2.setId("numOfIt");
         L2.setText("Number Of Iterations:");
         L2.setLayoutX(8);
         L2.setLayoutY(100);
         L2.setPrefHeight(28);
         L2.setPrefWidth(150);
         L2.setStyle(" -fx-font-weight: bold; -fx-background-radius:5;");
         L2.setAlignment(Pos.CENTER);
         
         T1.setId("v1");
         T1.setLayoutX(154);
         T1.setLayoutY(100);
         T1.setPrefHeight(28);
         T1.setPrefWidth(45);
         T1.setStyle("-fx-background-color:  #33887799;");
         //creating relative error text input 
         L3.setId("error");
         L3.setText("Absolute Relative Error:");
         L3.setLayoutX(214);
         L3.setLayoutY(100);
         L3.setPrefHeight(28);
         L3.setPrefWidth(200);
         L3.setStyle(" -fx-font-weight: bold; -fx-background-radius:5;");
         L3.setAlignment(Pos.CENTER);
         
         T2.setId("v2");
         T2.setLayoutX(389);
         T2.setLayoutY(100);
         T2.setPrefHeight(28);
         T2.setPrefWidth(60);
         T2.setStyle("-fx-background-color: #33887799;");
         
         ap.getChildren().addAll(L1,L2,L3,T1,T2);
}
	public void solve2(int n) {
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if(node.getId().equals("xs") || node.getId().equals("rt") || node.getId().equals("reachedIt"))
				return true;
			return false;
		});
		int m = 1;
		TextField tfError = (TextField) ap.lookup("#v2");
		if(!tfError.getText().equals(""))
			er = Double.parseDouble(tfError.getText());
		if(!signDigits.getText().equals(""))
			significantDigits = Integer.parseInt(signDigits.getText());
		if(method.equals("MNR2")) {
		TextField multiplicity = (TextField) ap.lookup("#mm");
		if(!multiplicity.getText().equals(""))
			m = Integer.parseInt(multiplicity.getText());
		else {
			showAlert("YOU MUST ENTER THE MULTIPLICITY");
			return;
		}
		}
		TextField tfMaxIt = (TextField) ap.lookup("#v1");
		if(!tfMaxIt.getText().equals(""))
			iter = Integer.parseInt(tfMaxIt.getText());
		else {
			showAlert("YOU MUST ENTER THE MAX NUMBER OF ITERATIONS");
			return;
		}
		Function<Double, Double> f = x -> 0.0;
		for(int i=0;i<n;i++) {
			String w = Integer.toString(i);
			String coeff = ((TextField) ap.lookup("#a"+w)).getText();
			String fun = "";
			if(((ComboBox)ap.lookup("#dl"+w)).getSelectionModel().getSelectedItem()!=null)
				fun = (String)((ComboBox)ap.lookup("#dl"+w)).getSelectionModel().getSelectedItem();
			else {
				showAlert("YOU MUST CHOOSE ALL FUNCTIONS");
				return;}
			String p = ((TextField) ap.lookup("#cc"+w)).getText();
			
			double x1 = 1.0,x2 = 1.0;
			Function<Double, Double> fToAdd = null;
			if(!coeff.equals(""))
				x1 = Double.parseDouble(coeff);
			if(!p.equals(""))
				x2 = Double.parseDouble(p);
			
			final double finalX1 = x1;
			final double finalX2 = x2;
			if(fun.equals("x^c")) {
				fToAdd = z -> finalX1*Math.pow(z,finalX2);
			}
			else if(fun.equals("e^(cx)")) {
				fToAdd = z -> finalX1*Math.exp(finalX2*z); 
			}
			else if(fun.equals("sin(cx)")) {
				fToAdd = z -> finalX1*Math.sin(finalX2*z);
			}
			else if(fun.equals("cos(cx)")) {
				fToAdd = z -> finalX1*Math.cos(finalX2*z);
			}
			Function<Double, Double> finalFToAdd = fToAdd;
			Function<Double, Double> tempF = f;
		    f = z -> tempF.apply(z) + finalFToAdd.apply(z);
		}
		
		double[] solution = new double[2];
		double lower;
		double upper;
		long startTime = System.nanoTime();
		if(method.equals("BIS") || method.equals("FP") || method.equals("SEC")) {
			TextField tff = (TextField) ap.lookup("#Lbound");
			TextField tff2 = (TextField) ap.lookup("#Ubound");
			if(tff.getText().equals("") || tff2.getText().equals("")){
				showAlert("YOU MUST ENTER INITIAL GUESSES");
				return;
			}
			lower = Double.parseDouble(tff.getText());
			upper = Double.parseDouble(tff2.getText());
			if(lower==upper) {
				showAlert("INITIAL BOUNDS CANNOT BE EQUAL !");
				return;
			}
			
			if(method.equals("BIS"))
				solution = Bisection.bisection(lower, upper,er , iter,significantDigits, f);
			else if(method.equals("FP")) 
				solution = FalsePosition.falsePosition(lower, upper,er , iter,significantDigits, f);
			else if(method.equals("SEC"))
				solution = secant.secant(f,lower,upper, iter, significantDigits, er);
		}
		else if(method.equals("FXP") || method.equals("NR") || method.equals("MNR")|| method.equals("MNR2")) {
			TextField tff = (TextField) ap.lookup("#Lbound");
			if(tff.getText().equals("")){
				showAlert("YOU MUST ENTER INITIAL GUESS");
				return;
			}
			lower = Double.parseDouble(tff.getText());
			if(method.equals("FXP"))
				solution = FixedPoint.fixedPoint(f, lower, iter, significantDigits, er);
			else if(method.equals("NR")) {
				solution = NewtonRaphason.NewtonRaphson(lower, er, iter, significantDigits, f);
			}
			else if(method.equals("MNR")) {
				solution = ModifiedNewtonRaphson.modifiedNewtonRaphson(lower, er, iter, significantDigits, f);
			}
			else if(method.equals("MNR2")) {
				solution = ModifiedNewtonRaphson.modifiedNewtonRaphson2(lower,m, er, iter, significantDigits, f);
			}
		}
		if(solution[0] == Double.POSITIVE_INFINITY) {
			showAlert("METHOD DIVERGES:(");
			return;
		}
		else if(Double.isNaN(solution[0])) {
			showAlert("NO ROOT FOUND :( \nTRY TO CHANGE THE METHOD OR INITIAL GUESS");
			return;
		}
		long endTime = System.nanoTime();
		double runTime = (endTime - startTime)/1000000.0 ; // runtime of the chosen method in milliseconds 
		 Label L = new Label();
         L.setId("xs");
         L.setText("Root = "+Double.toString(solution[0]));
         L.setLayoutX(35);
         L.setLayoutY(320);
         L.setPrefHeight(28);
         L.setPrefWidth(significantDigits*7.5+45);
         L.setStyle("-fx-background-color:  #33887799;");
         L.setAlignment(Pos.CENTER);
         ap.getChildren().add(L);
         Label LL = new Label();
         LL.setId("rt");
         LL.setText("Run Time = "+Double.toString(runTime)+"ms");
         LL.setLayoutX(35);
         LL.setLayoutY(350);
         LL.setPrefHeight(28);
         LL.setPrefWidth(250);
         LL.setStyle("-fx-background-color:  #33887799;");
         LL.setAlignment(Pos.CENTER);
         ap.getChildren().add(LL);
         Label LL2 = new Label();
         LL2.setId("reachedIt");
         LL2.setText("Reached Iterations = "+Integer.toString((int)solution[1]));
         LL2.setLayoutX(35);
         LL2.setLayoutY(380);
         LL2.setPrefHeight(28);
         LL2.setPrefWidth(250);
         LL2.setStyle("-fx-background-color:  #33887799;");
         LL2.setAlignment(Pos.CENTER);
         ap.getChildren().add(LL2);
	}
	public void solve(){
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if(node.getId().equals("xs") || node.getId().equals("rt") || node.getId().equals("reachedIt"))
				return true;
			return false;
		});
		
		double[][] coeff = new double[N][N];
		double[] consts = new double[N];
		double[] guess = new double[N];
		
		//initializing the coefficients matrix A
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				String w = Integer.toString(i);
	            String q = Integer.toString(j);
	            String lookId = "#a" + w + q;
	            TextField tf = (TextField) ap.lookup(lookId);
				if(tf.getText()=="")
					coeff[i][j] = 0;
	            else
	            	coeff[i][j] = Double.parseDouble(tf.getText());
			}
		}
		//initializing the constants matrix B 
		for(int i=0;i<N;i++) {
			String w = Integer.toString(i);
            String lookId = "#b" + w;
            TextField tf = (TextField) ap.lookup(lookId);
            if(tf.getText()=="")
            	consts[i] = 0;
            else
            	consts[i] = Double.parseDouble(tf.getText());
		}
		if(method.equals("")) {
			showAlert("NO METHOD IS SELECTED, PLEASE SELECT ONE.");
			return;
		}
		// initializing the initial guess vector ,max number of iterations and relative error
		else if(method.equals("JI") || method.equals("GS")) {
			for(int i=0;i<N;i++) {
				String w = Integer.toString(i);
				String lookId = "#ia" + w;
	            TextField tf = (TextField) ap.lookup(lookId);
	            if(tf.getText()=="")
	            	guess[i] = 0;
	            else
	            	guess[i] = Double.parseDouble(tf.getText());
			}
			TextField tf1 = (TextField) ap.lookup("#v1");
			TextField tf2 = (TextField) ap.lookup("#v2");
			if(!tf1.getText().equals(""))
				iter = Integer.parseInt(tf1.getText());
			if(!tf2.getText().equals(""))
				er = Double.parseDouble(tf2.getText());
		}
		// checking the System
		if(Solver.checkSolutionType(coeff, consts)==0) {
			showAlert("SYSTEM HAS NO SOLUTION :(");;
			return;
		}
		else if(Solver.checkSolutionType(coeff, consts)==2) {
			showAlert("SYSTEM HAS INFINITE SOLUTIONS..!!!");
			return;
		}
		
		double[] solution = new double[N];
		long startTime = System.nanoTime();
		long endTime = 0;
	    if(method.equals("GE"))
			solution = GaussElimination.gauss(N,coeff, consts, significantDigits);
		else if(method.equals("GJ"))
			solution = GuassJordanElimination.gaussJordan(N,coeff, consts, significantDigits);
		else if(method.equals("LUdo"))
		    solution = LUDoolittle.solveLUdo(coeff, consts, significantDigits);
		else if(method.equals("LUcr"))
			solution = LUCrout.solveLUcr(coeff, consts, significantDigits);
		else if(method.equals("LUch")) {
			if(Solver.forSymmetry(coeff))
				solution = LUCholesky.LUch(coeff, consts,significantDigits);
			else {
				showAlert("COEFFICIENTS MATRIX MUST BE SYMMETRIC !");
				return;
			}
		}
		else if(method.equals("JI")) {
			solution = Jacobi.jacobi(coeff, consts,guess,er,iter,significantDigits);
			endTime = System.nanoTime();
			if(!Solver.forDiagonalDominance(coeff)) {
				showAlert("THIS SYSTEM IS NOT DIAGONALLY DOMINANT,SO SOLUTION WILL NOT CONVERGE.");
			}
		}
		else if(method.equals("GS")) {
			solution = GaussSeidel.gSediel(coeff, consts,guess,er,iter,significantDigits);
			endTime = System.nanoTime();
			if(!Solver.forDiagonalDominance(coeff)) {
				showAlert("THIS SYSTEM IS NOT DIAGONALLY DOMINANT,SO SOLUTION MAY NOT CONVERGE.");
			}
		}
	    if(!method.equals("LUch") && !method.equals("JI") &&!method.equals("GS"))
	    	endTime = System.nanoTime();			
		double runTime = (endTime - startTime)/1000000.0 ; // runtime of the chosen method in milliseconds 
			
		int x = 35;
		double width =significantDigits*7.5+45;
		for(int i=0;i<N;i++){
				 Label L = new Label();
		         L.setId("xs");
		         L.setText("X"+Integer.toString(i)+" = "+Double.toString(solution[i]));
		         L.setLayoutX(x);
		         L.setLayoutY((N+1)*28 + 200 + 10*(N+1));
		         L.setPrefHeight(28);
		         L.setPrefWidth(width);
		         L.setStyle("-fx-background-color:  #33887799;");
		         L.setAlignment(Pos.CENTER);
		         ap.getChildren().add(L);
		         x += width+10;
		}
		Label L = new Label();
        L.setId("rt");
        L.setText("Run Time = "+Double.toString(runTime)+"ms");
        L.setLayoutX(35);
        L.setLayoutY((N+1)*28 + 250 + 10*(N+1));
        L.setPrefHeight(28);
        L.setPrefWidth(250);
        L.setStyle("-fx-background-color:  #33887799;");
        L.setAlignment(Pos.CENTER);
        ap.getChildren().add(L);
	}
	@FXML
	//selecting which method to solve with
	public void select() {
		ap.getChildren().removeIf(node -> {
			if(node.getId()==null)
				return false;
			else if(node.getId().contains("ix") || 
					node.getId().contains("ia")||
					node.getId().equals("numOfIt")||
					node.getId().equals("stopCondition")||
					node.getId().equals("error")||
					node.getId().equals("v1")||
					node.getId().equals("v2")
					)
				return true;
			return false;
		});
		String s = comBox.getSelectionModel().getSelectedItem();
		if(s.equals("Gauss Elimination"))
			method = "GE";
		else if(s.equals("Gauss-Jordan"))
			method = "GJ";
		else if(s.equals("LU Decomposition(Crout)"))
			method = "LUcr";
		else if(s.equals("LU Decomposition(Cholesky)"))
			method = "LUch";
		else if(s.equals("LU Decomposition(Doolittle)"))
			method = "LUdo";
		else if(s.equals("Gauss-Seidel"))
			method = "GS";
		else if(s.equals("Jacobi-Iteration"))
			method = "JI";
		else if(s.equals("Bisection"))
			method = "BIS";
		else if(s.equals("False-Position"))
			method = "FP";
		else if(s.equals("Fixed Point"))
			method = "FXP";
		else if(s.equals("Secant"))
			method = "SEC";
		else if(s.equals("Newton Raphason"))
			method = "NR";
		else if(s.equals("Modified Newton Raphason 1"))
			method = "MNR";
		else if(s.equals("Modified Newton Raphason 2"))
			method = "MNR2";
		
		if(!method.equals("JI") && !method.equals("GS")) {
			ap.getChildren().removeIf(node -> {
				if(node.getId()==null)
					return false;
				else if(node.getId().contains("ix") || node.getId().contains("ia"))
					return true;
				return false;
			});
		}
		else if(method.equals("JI") || method.equals("GS") && !tfNumOfEq.getText().equals(""))
			create2(Integer.parseInt(tfNumOfEq.getText()));
		
		if(method.equals("FP") || method.equals("BIS") ||  
				   method.equals("FXP") ||method.equals("SEC")||
				   method.equals("NR")||method.equals("MNR")||method.equals("MNR2")) {
					create4();
					}
	}
	private void showAlert(String s) {
        Alert alert = new Alert(AlertType.ERROR);
        alert.setTitle("Error Message");
        alert.setHeaderText(null);
        alert.setContentText(s);
        alert.showAndWait();// Show the alert
    }
}
