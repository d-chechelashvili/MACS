//
//  InputCell.swift
//  Login Screen
//
//  Created by dato che on 11/22/20.
//

import UIKit

class InputCell: UITableViewCell {
    
    @IBOutlet var inputField: UITextField!
    @IBOutlet var headerLabel: UILabel!
    @IBOutlet var visualHeader: UIImageView!
    @IBOutlet var readableButton: UIButton!
    @IBOutlet var wrapperView: UIView!

    override func awakeFromNib() {
        super.awakeFromNib()
        wrapperView.layer.cornerRadius = wrapperView.bounds.height/5
        wrapperView.layer.borderWidth = wrapperView.bounds.height/80
        wrapperView.layer.borderColor = UIColor.lightGray.cgColor
        inputField.borderStyle = .none
        // Initialization code
    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)
        // Configure the view for the selected state
    }
    
    @IBAction func readabilityChanged(){
        if inputField.isSecureTextEntry{
            readableButton.setBackgroundImage(UIImage(systemName: "eye.fill"), for: .normal)
        }else {
            readableButton.setBackgroundImage(UIImage(systemName: "eye.slash.fill"), for: .normal)
        }
        inputField.isSecureTextEntry = !inputField.isSecureTextEntry
    }
    
}
