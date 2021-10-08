//
//  ButtonCell.swift
//  Burger Menu
//
//  Created by dato che on 12/8/20.
//

import UIKit

class ButtonCellModel{
    var text: String
    
    init(text: String) {
        self.text = text
    }
}


class ButtonCell: UITableViewCell {

    @IBOutlet var buttonCellLabel: UILabel!
    
    private var model: ButtonCellModel!
    
    override func awakeFromNib() {
        super.awakeFromNib()
        // Initialization code
    }
    
    func configure(with model: ButtonCellModel){
        buttonCellLabel.text = model.text

    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)

        // Configure the view for the selected state
    }
    
    @IBAction func buttontapped(){
        buttonCellLabel.alpha = 0.3
        //buttonCellLabel.tintColor = .orange
        buttonCellLabel.textColor = .orange  // ver vxvdebi ratom ar cvlis fers
        UIView.animate(
            withDuration: 0.8,
            delay: 0.0,
            options: [.curveLinear],
            animations: {
                self.buttonCellLabel.alpha = 1.0
                self.buttonCellLabel.textColor = .black
                //self.buttonCellLabel.tintColor = .black
            },
            completion: nil
        )
    }
    
}

