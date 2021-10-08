//
//  ButtonCell.swift
//  Login Screen
//
//  Created by dato che on 11/22/20.
//

import UIKit

class ButtonCell: UITableViewCell {

    @IBOutlet var button: UIButton!
    var clickHandler: (() -> ())?
    
    override func awakeFromNib() {
        super.awakeFromNib()
        button.layer.cornerRadius = button.bounds.height/2.5
        button.layer.borderWidth = button.bounds.height/30
        button.layer.borderColor = UIColor.blue.cgColor
        button.titleLabel?.font = .systemFont(ofSize: 20)
        // Initialization code
    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)

        // Configure the view for the selected state
    }
    
    @IBAction func buttonClicked(){
        clickHandler?()
    }
    
}
