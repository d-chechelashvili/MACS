//
//  TableController.swift
//  Login Screen
//
//  Created by dato che on 11/22/20.
//

import UIKit
//import WXWaveView

class TableController: UIViewController{
    
   @IBOutlet var mainTableView: UITableView!
    
    override func viewDidLoad(){
        super.viewDidLoad()
        listenToKeyboardNotifications()
        mainTableView.delegate = self
        mainTableView.dataSource = self
        mainTableView.register(UINib(nibName: "ButtonCell", bundle: nil),
                               forCellReuseIdentifier: "ButtonCell")
        mainTableView.register(UINib(nibName: "ForgotPasswordCell", bundle: nil),
                               forCellReuseIdentifier: "ForgotPasswordCell")
        mainTableView.register(UINib(nibName: "InputCell", bundle: nil),
                               forCellReuseIdentifier: "InputCell")
        mainTableView.register(UINib(nibName: "LoginViewCell", bundle: nil),
                               forCellReuseIdentifier: "LoginViewCell")
        
        mainTableView.separatorStyle = .none
        mainTableView.allowsSelection = false
        mainTableView.keyboardDismissMode = .interactive
    }
    
    func listenToKeyboardNotifications(){
        NotificationCenter.default.addObserver(
            self,
            selector: #selector(handleKeyboardShow(notification:)),
            name: UIApplication.keyboardWillShowNotification,
            object: nil)
    
        NotificationCenter.default.addObserver(
            self,
            selector: #selector(handleKeyboardHide(notification:)),
            name: UIApplication.keyboardWillHideNotification,
            object: nil)
    }
    
    @objc func handleKeyboardShow(notification: Notification){
        guard let userInfo = notification.userInfo,
              let duration = userInfo[UIResponder.keyboardAnimationDurationUserInfoKey] as? Double,
              let endFarme = userInfo[UIResponder.keyboardFrameEndUserInfoKey] as? NSValue
        else {return}
        
        UIView.animate(
            withDuration: duration,
            animations: {
                self.mainTableView.contentInset = UIEdgeInsets(
                    top: 0,
                    left: 0,
                    bottom: endFarme.cgRectValue.height - self.view.safeAreaInsets.bottom,
                    right: 0)
            })
        
        self.mainTableView.scrollToRow(
            at: IndexPath(row: 1, section: 3),
            at: .bottom,
            animated: false
        )
    }
    
    @objc func handleKeyboardHide(notification: Notification){
        guard let userInfo = notification.userInfo,
              let duration = userInfo[UIResponder.keyboardAnimationDurationUserInfoKey] as? Double
        else {return}
        
        
        UIView.animate(
            withDuration: duration,
            animations: {
                self.mainTableView.contentInset = UIEdgeInsets(
                    top: 0,
                    left: 0,
                    bottom: 0,
                    right: 0)
            })
    }
}

extension TableController: UITableViewDelegate, UITableViewDataSource{
    func numberOfSections(in tableView: UITableView) -> Int {
        return 4
    }
    
    func tableView(_ tableView: UITableView, numberOfRowsInSection section: Int) -> Int {
        if section == 0 || section == 2 {
            return 1
        }
        return 2
    }
    
    func tableView(_ tableView: UITableView, heightForRowAt indexPath: IndexPath) -> CGFloat {
        if indexPath.section == 0 {
            return UIScreen.main.bounds.height * 0.4
        }else if indexPath.section == 1 {
            return UIScreen.main.bounds.height * 0.11
        }else if indexPath.section == 2 {
            return 44
        }
        return UIScreen.main.bounds.height * 0.12
    }
    
    func tableView(_ tableView: UITableView, cellForRowAt indexPath: IndexPath) -> UITableViewCell {
        if indexPath.section == 0 {
            return initLoginViewCell(tableView, indexPath)
        }else if indexPath.section == 1 {
            return initInputCell(tableView, indexPath)
        }else if indexPath.section == 2 {
           return initForgotPasswordCell(tableView, indexPath)
        }
        return initButtonCell(tableView, indexPath)
    }
    
    func initLoginViewCell(_ tableView: UITableView,_ indexPath: IndexPath) -> UITableViewCell{
        let cell = mainTableView.dequeueReusableCell(withIdentifier: "LoginViewCell", for: indexPath)
        
        //არ მუშაობს რატომღაც, 3 საათია ვეწვალები, მოვიდა დედლაინი და ვაკომენტარებ, არ გამომივიდა :(
        //credits: https://github.com/WelkinXie/WXWaveView

//        if let loginViewCell = cell as? LoginViewCell{
//        let waveView = WXWaveView.add(to: loginViewCell.headerView,
//                                      withFrame: CGRect(x: 0,
//                                                        y: loginViewCell.headerView.bounds.height -10,
//                                                        width: loginViewCell.headerView.bounds.width,
//                                                        height: 10))
//            waveView?.waveColor = .white
//            waveView?.wave()
//        }

        return cell
    }
    
    func initInputCell(_ tableView: UITableView,_ indexPath: IndexPath) -> UITableViewCell{
        let cell = mainTableView.dequeueReusableCell(withIdentifier: "InputCell", for: indexPath)
        if let inputCell = cell as? InputCell{
            if indexPath.row == 0 {
                inputCell.visualHeader.image = UIImage(systemName: "envelope")
                inputCell.headerLabel.text = "Email"
                inputCell.inputField.placeholder = "test@gmail.com"
                inputCell.readableButton.isHidden = true
            }else{
                inputCell.visualHeader.image = UIImage(systemName: "lock")
                inputCell.inputField.isSecureTextEntry = true
                inputCell.headerLabel.text = "Password"
                inputCell.inputField.placeholder = "Hello!"
            }
        }
        return cell
    }
    
    func initForgotPasswordCell(_ tableView: UITableView,_ indexPath: IndexPath) -> UITableViewCell{
        let cell = mainTableView.dequeueReusableCell(withIdentifier: "ForgotPasswordCell", for: indexPath)
        return cell
    }
    
    func initButtonCell(_ tableView: UITableView,_ indexPath: IndexPath) -> UITableViewCell{
        let cell = mainTableView.dequeueReusableCell(withIdentifier: "ButtonCell", for: indexPath)
        if let buttonCell = cell as? ButtonCell{
            if indexPath.row == 0 {
                initLoginCell(buttonCell: buttonCell)
            }else{
                buttonCell.button.setTitle("Sign Up", for: .normal)
                buttonCell.button.setTitleColor(.blue, for: .normal)
            }
        }
        return cell
    }
    
    func initLoginCell(buttonCell: ButtonCell){
        buttonCell.clickHandler = {
            let cell1 = self.mainTableView.cellForRow(at: IndexPath(row: 0, section: 1))
            let cell2 = self.mainTableView.cellForRow(at: IndexPath(row: 1, section: 1))
            
            guard let emailCell = cell1 as? InputCell,
                  let passwordCell = cell2 as? InputCell
            else {return}
                
            if emailCell.inputField.text == "test@mail.com" &&
                passwordCell.inputField.text == "1234" {
                
                let mainStoryBoard = UIStoryboard(name: "Main", bundle: nil)
                let dialPadTBC = mainStoryBoard.instantiateViewController(withIdentifier: "DialPadTBC")
                dialPadTBC.modalPresentationStyle = .overFullScreen
                self.present(dialPadTBC, animated: true, completion: nil)
                //print("CHEMI GAMOQVABULIA")
            }else{
                //print("SXVAGAN MOVXVDI")
            }
        }
        
        buttonCell.button.setTitle("Login", for: .normal)
        buttonCell.button.setTitleColor(.white, for: .normal)
        buttonCell.button.backgroundColor = .blue
    }
    
}
