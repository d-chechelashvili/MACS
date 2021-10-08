//
//  ForecastViewController.swift
//  Weather App
//
//  Created by dato che on 2/7/21.
//

import UIKit

class ForecastSection{
    var forecast = [TimeStamp]()
    
    var numberOfRows: Int{
        return forecast.count
    }
    
    init(forecast: [TimeStamp]) {
        self.forecast = forecast
    }
}

class ForecastViewController: UIViewController{
    @IBOutlet var tableView: UITableView!
    @IBOutlet var activiyIndicator: UIActivityIndicatorView!
    @IBOutlet var activityIndicatorPage: UIView!
    @IBOutlet var errorPage: UIView!
    @IBOutlet var reloadButton: UIButton!
    
    private let service = Service(for: .prolonged)
    
    private var tableData = [ForecastSection]()
    
     
    override func viewDidLayoutSubviews() {
        super.viewDidLayoutSubviews()
        reloadButton.layer.cornerRadius = reloadButton.frame.width / 10
    }
    
    override func viewDidLoad(){
        super.viewDidLoad()
        makeBarsTransparent()
        setupTableView()
     }
    
    override func viewDidAppear(_ animated: Bool) {
        loadForecast()
    }
    
    @IBAction func handleRefresh(){
        loadForecast()
    }
    
    private func makeBarsTransparent() {
        self.navigationController?.navigationBar.setBackgroundImage(UIImage(), for: .default)
        self.navigationController?.navigationBar.shadowImage = UIImage()
        self.navigationController?.navigationBar.layoutIfNeeded()
    }
    
    
    private func setupTableView(){
        tableView.register(UINib(nibName: "ForecastCell", bundle: nil), forCellReuseIdentifier: "ForecastCell")
        
        tableView.register(UINib(nibName: "ForecastHeader", bundle: nil), forHeaderFooterViewReuseIdentifier: "ForecastHeader")
        
        tableView.delegate = self
        tableView.dataSource = self
        tableView.separatorStyle = .none
        tableView.allowsSelection = false
        tableView.backgroundView?.backgroundColor = .clear
        
    }
    
    private func loadForecast(){
        tableData.removeAll()
        errorPage.isHidden = true
        guard let cityName = getCityName() else {return}
        activityIndicatorPage.isHidden = false
        service.loadForecast(for: cityName) { [weak self] result in
            guard let self = self else {return}
            DispatchQueue.main.async {
                switch result{
                case .success(let forecast):
                    let timestamps = forecast.list
                    self.parseForecast(timestamps)
                    self.activityIndicatorPage.isHidden = true
                    self.tableView.reloadData()
                case .failure(_):
                    self.activityIndicatorPage.isHidden = true
                    self.errorPage.isHidden = false
                }
            }
        }
    }
    
    private func getCityName() -> String?{
        if let viewControllers = self.tabBarController?.viewControllers{
            if let todayNC = viewControllers[0] as? UINavigationController{
                if let todayVC = todayNC.viewControllers[0] as? TodayViewController{
                    return todayVC.getSelectedCity()
                }
            }
        }
        return nil
    }
    
    private func parseForecast(_ timestamps: [TimeStamp]){
        var i = 0
        var currDay = timestamps[i].dt_txt.dropLast(9)
        var currForecast: [TimeStamp] = []
        while i < timestamps.count{
            if !timestamps[i].dt_txt.hasPrefix(currDay){
                currDay = timestamps[i].dt_txt.dropLast(9)
                let todayForecast = ForecastSection(forecast: currForecast)
                tableData.append(todayForecast)
                currForecast.removeAll()
            }
            currForecast.append(timestamps[i])
            i += 1
        }
    }
}

extension ForecastViewController: UITableViewDelegate, UITableViewDataSource{
    
    func numberOfSections(in tableView: UITableView) -> Int {
        return tableData.count
    }
    func tableView(_ tableView: UITableView, numberOfRowsInSection section: Int) -> Int {
        return tableData[section].numberOfRows
    }
    
    func tableView(_ tableView: UITableView, cellForRowAt indexPath: IndexPath) -> UITableViewCell {

        let cell = tableView.dequeueReusableCell(withIdentifier: "ForecastCell", for: indexPath)
        if let forecastCell = cell as? ForecastCell{
            let timestamp = tableData[indexPath.section].forecast[indexPath.row]
            forecastCell.descriptionLabel.text = timestamp.weather[0].description
            forecastCell.temperatureLabel.text = Int(round(timestamp.main.temp - 273.15)).description + "°C"
            forecastCell.timeLabel.text = timestamp.dt_txt.dropFirst(11).description.dropLast(3).description
            forecastCell.weatherImage.image = getWeatherImage(timestamp.weather[0].icon)
        }
        return cell
    }
    
    func getWeatherImage(_ icon: String) -> UIImage{
        let imageURl = "https://openweathermap.org/img/w/" + icon + ".png"
        guard let url = URL(string: imageURl) else { return UIImage() }
        if let data = try? Data(contentsOf: url){
            if let image = UIImage(data: data){
                return image
            }
        }
        return UIImage()
    }
    
    
    func tableView(_ tableView: UITableView, viewForHeaderInSection section: Int) -> UIView? {
        let header = tableView.dequeueReusableHeaderFooterView(withIdentifier: "ForecastHeader")
        let timestamp = tableData[section].forecast[0]
        let currDate = String(timestamp.dt_txt.dropLast(9))
        
        let dateFormatter = DateFormatter()
        dateFormatter.dateFormat = "yyyy-MM-dd"
        let date = dateFormatter.date(from: currDate)
        dateFormatter.dateFormat = "EEEE"
        if let forecastHeader = header as? ForecastHeader{
            forecastHeader.titleLabel.text = dateFormatter.string(from: date!)
        }
        
        return header
    }
    
    func tableView(_ tableView: UITableView, heightForHeaderInSection section: Int) -> CGFloat {
        return 44
    }
    
    func tableView(_ tableView: UITableView, heightForRowAt indexPath: IndexPath) -> CGFloat {
        return 72
    }
}


