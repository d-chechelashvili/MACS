//
//  ViewController.swift
//  Dial Pad
//
//  Created by dato che on 11/11/20.
//

import UIKit

class ViewController: UIViewController {
    
    @IBOutlet var number: UILabel!
    @IBOutlet var dialPadRow1: DialPadRow!
    @IBOutlet var dialPadRow2: DialPadRow!
    @IBOutlet var dialPadRow3: DialPadRow!
    @IBOutlet var dialPadRow4: DialPadRow!
    
    override func viewDidLoad() {
        super.viewDidLoad()
        setupRows()
        // Do any additional setup after loading the view.
    }
    
    override func viewDidLayoutSubviews() {
        super.viewDidLayoutSubviews()
        setupButtons()
    }

    
    func setupRows(){
        dialPadRow1.delegate = self
        dialPadRow1.setup()
        dialPadRow2.delegate = self
        dialPadRow2.setup()
        dialPadRow3.delegate = self
        dialPadRow3.setup()
        dialPadRow4.delegate = self
        dialPadRow4.setup()
    }
    
    func setupButtons(){
        setupButton(button: dialPadRow1.dialPadButton1, digit: "1", text: "")
        setupButton(button: dialPadRow1.dialPadButton2, digit: "2", text: "A B C")
        setupButton(button: dialPadRow1.dialPadButton3, digit: "3", text: "D E F")
        setupButton(button: dialPadRow2.dialPadButton1, digit: "4", text: "G H I")
        setupButton(button: dialPadRow2.dialPadButton2, digit: "5", text: "J K L")
        setupButton(button: dialPadRow2.dialPadButton3, digit: "6", text: "M N O")
        setupButton(button: dialPadRow3.dialPadButton1, digit: "7", text: "P Q R S")
        setupButton(button: dialPadRow3.dialPadButton2, digit: "8", text: "T U V")
        setupButton(button: dialPadRow3.dialPadButton3, digit: "9", text: "W X Y Z")
        setupButton(button: dialPadRow4.dialPadButton1, digit: "*", text: "")
        setupButton(button: dialPadRow4.dialPadButton2, digit: "0", text: "+")
        setupButton(button: dialPadRow4.dialPadButton3, digit: "#", text: "")
    }
    
    func setupButton(button: DialPadButton, digit: String, text: String){
        button.digitLabel.text = digit
        button.textLabel.text = text
        button.setup()
    }

    
    @IBAction func deleteButtonTapped(){
        if(number.text != nil && number.text != ""){
            number.text = String(number.text!.dropLast())
        }
    }
    
    @IBAction func addNumberButtonTapped(){
        let mainStoryboard = UIStoryboard(name: "Main", bundle: nil)

        if let addNumberVC = mainStoryboard.instantiateViewController(withIdentifier: "AddNumberVC") as? UIViewController {
           navigationController?.pushViewController(addNumberVC, animated: true)
        }
    }
    
    @IBAction func dialButtonTapped(){
        let mainStoryboard = UIStoryboard(name: "Main", bundle: nil)
        
        if let dialVCNavigator = mainStoryboard.instantiateViewController(withIdentifier: "DialVCNavigation") as? DialVC {
            
            dialVCNavigator.modalPresentationStyle = .overFullScreen
            present(dialVCNavigator, animated: true, completion: nil)
        }
    }
}

extension ViewController: DialPadButtonDelegate{
    func dialPadButtonTapped(_ sender: DialPadButton){
        number.text?.append(sender.digitLabel.text!)
    }
}

