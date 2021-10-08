//
//  DialPadButton.swift
//  Dial Pad
//
//  Created by dato che on 11/11/20.
//

import 	UIKit

protocol DialPadButtonDelegate: AnyObject{
    func dialPadButtonTapped(_ sender:DialPadButton)
}

class DialPadButton: BaseReusableView{
    @IBOutlet var digitLabel: UILabel!
    @IBOutlet var textLabel: UILabel!
    @IBOutlet var dialPadButton: UIButton!
    
    
    weak var delegate: DialPadButtonDelegate?
    
    override func setup() {
        contentView.layer.cornerRadius = contentView.frame.width/2
    }
    
    
    
    @IBAction func handleDialPadButtonTap() {
        delegate?.dialPadButtonTapped(self)
    }
    
}
