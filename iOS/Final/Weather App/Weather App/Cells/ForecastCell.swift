//
//  ForecastCell.swift
//  Weather App
//
//  Created by dato che on 2/7/21.
//

import UIKit

class ForecastCell: UITableViewCell {

    @IBOutlet var weatherImage: UIImageView!
    @IBOutlet var timeLabel: UILabel!
    @IBOutlet var descriptionLabel: UILabel!
    @IBOutlet var temperatureLabel: UILabel!
        
    override func awakeFromNib() {
        super.awakeFromNib()
        // Initialization code
    }

    override func setSelected(_ selected: Bool, animated: Bool) {
        super.setSelected(selected, animated: animated)

        // Configure the view for the selected state
    }
    
}
