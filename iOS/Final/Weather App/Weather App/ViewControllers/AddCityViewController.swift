//
//  AddCityViewController.swift
//  Weather App
//
//  Created by dato che on 2/4/21.
//

import UIKit


protocol AddCityViewControllerDelegate: AnyObject{
    func addCityButtonTapped(_ sender:AddCityViewController, _ serviceResponse: ServiceResponse)
}

class AddCityViewController: UIViewController{
    
    @IBOutlet var cityNameField: UITextField!
    @IBOutlet var addCityView: UIView!
    @IBOutlet var errorView: UIView!
    @IBOutlet var errorLabel: UILabel!
    @IBOutlet var addCityButton: UIButton!
    @IBOutlet var activityIndicator: UIActivityIndicatorView!
    @IBOutlet var activityIndicatorOuterView: UIView!
    @IBOutlet var yPosConstraint: NSLayoutConstraint!
    
    weak var delegate: AddCityViewControllerDelegate?
    
    private var service = Service(for: .today)
    
    override func viewDidLayoutSubviews() {
        super.viewDidLayoutSubviews()
        addCityView.layer.cornerRadius = addCityView.frame.height / 10
        errorView.layer.cornerRadius = errorView.frame.width / 15
        activityIndicatorOuterView.layer.cornerRadius = activityIndicatorOuterView.frame.width/2
    }
    
    override func viewDidLoad() {
        super.viewDidLoad()
        listenToKeyboardNotifications()
    }
    
    func listenToKeyboardNotifications(){
        cityNameField.delegate = self
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
                self.yPosConstraint.constant = endFarme.cgRectValue.height - self.view.safeAreaInsets.bottom + 10
                self.yPosConstraint.priority = UILayoutPriority(.init(1000))
           })
    }
    
    @objc func handleKeyboardHide(notification: Notification){
        guard let userInfo = notification.userInfo,
              let duration = userInfo[UIResponder.keyboardAnimationDurationUserInfoKey] as? Double
        else {return}
        
        
        UIView.animate(
            withDuration: duration,
            animations: {
                self.yPosConstraint.priority = UILayoutPriority(.init(500))
            })
    }
    
    private var group = DispatchGroup()
    @IBAction func handleAddCity(){
        self.activityIndicatorOuterView.isHidden = false
        service.loadWeather(for: cityNameField.text ?? ""){ [weak self] result in
            guard let self = self else {return}
            DispatchQueue.main.async {
                switch result {
                case .success(let weather):
                    self.delegate?.addCityButtonTapped(self, weather)
                    self.dismiss(animated: true, completion: nil)
                case .failure(let error):
                    self.errorLabel.text = error.localizedDescription
                    self.errorView.isHidden = false
                }
                self.activityIndicatorOuterView.isHidden = true
            }
        }
    }
    
    @IBAction func handleOuterClick(){
        dismiss(animated: true, completion: nil)
    }
}

extension AddCityViewController: UITextFieldDelegate{
    func textFieldShouldReturn(_ textField: UITextField) -> Bool {
        textField.resignFirstResponder()
        return true
    }
}
